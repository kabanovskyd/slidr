# written by Matthew Shebet at the Broad Institute

# activate the project environment so all required packages are available
using Pkg
project_path = ENV["JULIA_PROJECT_PATH"]
Pkg.activate(project_path)
Pkg.precompile()

using CSV
using HDF5
using FASTX
using CodecZlib
using IterTools: product
using StatsBase: countmap, sample
using DataFrames
using StringViews
using LinearAlgebra: dot
using Combinatorics: combinations
using Dates


# Progress-logging vocabulary, deliberately identical to the Python side's (helpers.log_ts /
# log_detail / log_section / format_duration). This script's stdout is captured into
# spatial_barcodes.log, and a reader moving between that file and runtime.log should not have to
# adjust to a different layout: one timestamp column, details hanging beneath, a rule per phase.
const TS_WIDTH = 10

logts(msg) = (println(Dates.format(now(), "HH:MM:SS") * "  " * msg); flush(stdout))
logdetail(msg) = (println(" "^TS_WIDTH * msg); flush(stdout))

function logsection(title, width=74)
    prefix = "── $title "
    println()
    println(prefix * "─"^max(0, width - length(prefix)))
    flush(stdout)
end

"""
Render an elapsed time compactly: '45s', '12m34s', '2h07m' -- the same shapes format_duration()
produces on the Python side.
"""
function fmtduration(seconds)
    total = round(Int, seconds)
    hours, rem = divrem(total, 3600)
    minutes, secs = divrem(rem, 60)
    hours > 0 && return "$(hours)h$(lpad(minutes, 2, '0'))m"
    minutes > 0 && return "$(minutes)m$(lpad(secs, 2, '0'))s"
    return "$(secs)s"
end

# Time a phase: logs a start line, runs `f`, then logs a tick with the elapsed time. Phase timings are
# what make a slow run diagnosable, and they are cheap to record here where the work happens.
function timed(label, f)
    logts(label)
    t0 = time()
    result = f()
    logts("✓ $label — $(fmtduration(time() - t0))")
    return result
end


# Format a fatal-error message followed by actionable troubleshooting bullets, matching the
# "<message> / Troubleshooting: / • ..." layout the Python side of the pipeline uses. This script
# runs as a subprocess whose captured output (spatial_barcodes.log) is the only record of why it
# stopped, so a failed assertion has to say what to do about it and not just what went wrong.
function trouble(msg, hints...)
    io = IOBuffer()
    print(io, msg)
    if !isempty(hints)
        print(io, "\nTroubleshooting:")
        for hint in hints
            print(io, "\n • ", hint)
        end
    end
    return String(take!(io))
end


# args[1] — path to a directory containing all FASTQ files for this sample (R1 and R2 pairs)
# args[2] — path to a directory containing the puck CSV file(s) for this sample
# args[3] - sample ID string
# both directories must contain files for a single run only — all files in each are processed together
if length(ARGS) != 3
    error(trouble(
        "Usage: julia spatial_count.jl fastqpath puckpath sample_id",
        "This script is normally invoked by the pipeline's spatial-count stage, not by hand",
        "Run it through the pipeline with `./slidr --bcl <BCL_ID> --spatial-count`",
        "Got $(length(ARGS)) argument(s): $(ARGS)",
    ))
end
fastqpath = ARGS[1]
puckpath = ARGS[2]
sample_id = ARGS[3]

logsection("spatial barcode counting")
logts("sample $(sample_id)")
logdetail("fastqs → $(fastqpath)")
logdetail("puck   → $(puckpath)")

@assert isdir(fastqpath) trouble(
    "FASTQ path not found: $(fastqpath)",
    "Run the demultiplexing stage first with --mkfastq (or the whole pipeline with --run-all)",
    "If the reads were demultiplexed elsewhere, pass that directory with --fastqs",
)
@assert !isempty(readdir(fastqpath)) trouble(
    "FASTQ path is empty: $(fastqpath)",
    "The demultiplexing stage produced no spatial-barcode FASTQs for this sample",
    "Check the `SB Index` and `SB Lane` metadata columns match the sequencing run",
    "Re-run the demultiplexing stage with --mkfastq --force",
)

@assert isdir(puckpath) trouble(
    "Puck path not found: $(puckpath)",
    "The pipeline copies the puck CSV here before running this script, so this means that step did not happen",
    "Check `paths.puck_path` in the config file, and the `Puck ID` metadata column",
)
@assert !isempty(readdir(puckpath)) trouble(
    "Puck path is empty: $(puckpath)",
    "No puck CSV was staged or generated for this sample",
    "Check the `Puck ID` metadata column names a <puck id>.csv present in `paths.puck_path`",
    "Or set `paths.raw_barcodes_path` so the puck can be generated from BeadBarcodes.txt/BeadLocations.txt",
)

##### load the FASTQ data ######################################################

# split all files in the FASTQ directory into R1 and R2 lists by filename token
fastqs = readdir(fastqpath, join=true)
R1s = filter(s -> occursin(sample_id, s) && occursin("_R1_", s), fastqs)
R2s = filter(s -> occursin(sample_id, s) && occursin("_R2_", s), fastqs)
logts("found $(length(R1s)) FASTQ pair(s)")
for (r1, r2) in zip(R1s, R2s)
    logdetail("$(basename(r1))  +  $(basename(r2))")
end

# require at least one pair and matching filenames between R1 and R2
@assert length(R1s) == length(R2s) > 0 trouble(
    "expected an equal, non-zero number of R1 and R2 FASTQs for sample '$(sample_id)', but found $(length(R1s)) R1(s) and $(length(R2s)) R2(s) in $(fastqpath)",
    "Files are selected by containing both the sample ID and an `_R1_`/`_R2_` token in their name",
    "Check both reads of every lane are present and gzipped (.fastq.gz)",
    "A truncated or aborted demultiplexing run is the usual cause; re-run with --mkfastq --force",
    "Files present: $(basename.(readdir(fastqpath)))",
)
@assert [replace(R1, "_R1_"=>"", count=1) for R1 in R1s] == [replace(R2, "_R2_"=>"", count=1) for R2 in R2s] trouble(
    "the R1 and R2 FASTQs in $(fastqpath) do not pair up by filename",
    "Each R1 must have an R2 whose name is identical apart from the `_R1_`/`_R2_` token",
    "R1s: $(basename.(R1s))",
    "R2s: $(basename.(R2s))",
    "A leftover file from an earlier run is the usual cause; clear the directory and re-run with --mkfastq --force",
)

# return the length of the first read in a gzipped FASTQ file
function fastq_seq_len(path)
    reader = path |> open |> GzipDecompressorStream |> FASTQ.Reader
    try
        return(reader |> first |> FASTQ.sequence |> length)
    finally
        close(reader)
    end
end

# all R1s must have the same read length, and all R2s must have the same read length
# (mixing read lengths within one lane is almost always a sample prep error)
@assert length(unique([fastq_seq_len(R1) for R1 in R1s])) == 1 trouble(
    "the R1 FASTQs have different sequence lengths: $(unique([fastq_seq_len(R1) for R1 in R1s]))",
    "All reads processed together must share one read structure, or barcodes are sliced at the wrong offsets",
    "This normally means FASTQs from two different sequencing runs ended up in the same directory",
    "Check the `Merge Spatial From BCL` metadata column, and that $(fastqpath) holds only this run's reads",
    "R1s: $(basename.(R1s))",
)
@assert length(unique([fastq_seq_len(R2) for R2 in R2s])) == 1 trouble(
    "the R2 FASTQs have different sequence lengths: $(unique([fastq_seq_len(R2) for R2 in R2s]))",
    "All reads processed together must share one read structure, or barcodes are sliced at the wrong offsets",
    "This normally means FASTQs from two different sequencing runs ended up in the same directory",
    "Check the `Merge Spatial From BCL` metadata column, and that $(fastqpath) holds only this run's reads",
    "R2s: $(basename.(R2s))",
)

# scan the first 100k records to determine which file carries the UP sequence —
# in standard Slide-Tag layout R2 contains the spatial barcode + UP site, but
# some demultiplexers swap R1/R2, so we detect this automatically
function learn_switch(R1, R2)
    # both reads are indexed up to base 26 below (regardless of which one eventually turns out
    # to be the spatial read), so both — not just one of the two — must be long enough for that;
    # an OR here would let a too-short read slip through and crash the loop below with a BoundsError
    @assert (fastq_seq_len(R1) >= 26 && fastq_seq_len(R2) >= 26) trouble(
        "both R1 and R2 must have a read length of at least 26 bases to scan for the UP site, but are $(fastq_seq_len(R1)) and $(fastq_seq_len(R2)) bases",
        "Bases 9-26 of each read are compared against the UP sequence to work out which read carries the spatial barcode",
        "A read this short means the sequencing run used the wrong read configuration for Slide-Tag",
        "Check you are pointing at the spatial-barcode library and not the gene-expression one (the `SB Index` metadata column)",
        "Files: $(basename(R1)) / $(basename(R2))",
    )
    i1 = R1 |> open |> GzipDecompressorStream |> FASTQ.Reader
    i2 = R2 |> open |> GzipDecompressorStream |> FASTQ.Reader
    UPseq = String31("TCTTCAGCGTTCCCGAGA")
    s1 = 0 ; s2 = 0
    try
        # manual paired iteration (rather than zip(i1, i2), which silently truncates to the
        # shorter reader) so a mismatched/corrupted mate file is surfaced as a warning instead
        # of silently shrinking the sample used for swap detection
        i = 0
        next1 = iterate(i1)
        next2 = iterate(i2)
        while next1 !== nothing && next2 !== nothing && i < 100000
            i += 1
            record1 = next1[1]; record2 = next2[1]
            next1 = iterate(i1, next1[2])
            next2 = iterate(i2, next2[2])
            s1 += FASTQ.sequence(record1)[9:26] == UPseq
            s2 += FASTQ.sequence(record2)[9:26] == UPseq
        end
        if i < 100000 && (next1 === nothing) != (next2 === nothing)
            @warn "R1/R2 have mismatched read counts within the first 100k records sampled for swap detection ($R1 / $R2); swap-detection statistics below are based on fewer than 100k paired records"
        end
    finally
        close(i1)
        close(i2)
    end
    logdetail("$(basename(R1)): R1 matched UP $(s1)x, R2 matched UP $(s2)x")
    # return true if R1 is the spatial read (switched), false if R2 is (standard)
    return(s2 >= s1 ? false : true)
end

# Switch R1 and R2 if needed
logts("detecting which read carries the spatial barcode")
res_list = [learn_switch(R1, R2) for (R1, R2) in zip(R1s, R2s)]
@assert all(res_list) || !any(res_list) trouble(
    "the R1/R2 read assignment is not consistent across the FASTQ pairs (per-pair swap detection: $(res_list))",
    "In some pairs R1 carries the spatial barcode and in others R2 does, so one swap decision cannot be applied to all of them",
    "This normally means FASTQs from two differently demultiplexed runs ended up in the same directory",
    "Check the `Merge Spatial From BCL` metadata column, and that $(fastqpath) holds only this run's reads",
    "Split the mismatched pairs into separate runs, or re-demultiplex them consistently",
)
switch = all(res_list)

if switch == true
    logdetail("R1/R2 are swapped in these files — corrected")
    temp = R1s
    R1s = R2s
    R2s = temp
end

# after the (possible) swap above, `R2s` is always the file used as the spatial-barcode read for
# the rest of the run (record[2] below is sliced 1:32 for sb1/UP/sb2) — confirm that read is
# actually long enough, tied to the real post-swap assignment rather than an either/or guard
@assert fastq_seq_len(R2s[1]) >= 32 trouble(
    "the spatial-barcode read (R2, after R1/R2 swap correction) must have a read length of at least 32 bases, but is $(fastq_seq_len(R2s[1]))",
    "Bases 1-32 hold the two spatial-barcode halves and the UP site between them, so a shorter read cannot be decoded",
    "Check the sequencing run used the read configuration Slide-Tag expects",
    "Check you are pointing at the spatial-barcode library and not the gene-expression one (the `SB Index` metadata column)",
    "File: $(basename(R2s[1]))",
)
# the cell-barcode read (R1, after swap) is sliced 1:16 (cb) and 17:28 (umi) in process_fastqs, so
# it must be >= 28bp; assert it upfront too, symmetric with the spatial-read guard above, rather
# than letting a short R1 throw a BoundsError deep inside the read loop on the first record
@assert fastq_seq_len(R1s[1]) >= 28 trouble(
    "the cell-barcode read (R1, after R1/R2 swap correction) must have a read length of at least 28 bases (16bp cell barcode + 12bp UMI), but is $(fastq_seq_len(R1s[1]))",
    "Bases 1-16 hold the cell barcode and 17-28 the UMI, so a shorter read cannot be decoded",
    "Check the sequencing run used the read configuration Slide-Tag expects",
    "File: $(basename(R1s[1]))",
)

##### load the puck data #######################################################

# read all puck CSV files; columns are (sb, x, y) with no header
pucks = filter(x -> endswith(x, ".csv"), readdir(puckpath, join=true)) ; println("Pucks: ", basename.(pucks))
puckdfs = DataFrame[]
for puck in pucks
    raw = CSV.read(puck, DataFrame, header=false, types=[String15, Float64, Float64])
    # validate the column count before destructuring into named (sb, x, y) fields — puck CSVs
    # have no header, so a malformed file (wrong number of columns) would otherwise fail with a
    # confusing error from rename!, or silently misalign columns if counts happened to coincide
    @assert ncol(raw) == 3 trouble(
        "Puck CSV $(puck) must have exactly 3 columns (spatial barcode, x, y) with no header; found $(ncol(raw))",
        "Check the file has no header row -- a header is read as data and shifts everything",
        "Check the delimiter is a comma, not a tab or semicolon",
        "Inspect the first lines: `head -3 $(puck)`",
        "Delete the file and re-run to regenerate it from `paths.raw_barcodes_path`",
    )
    rename!(raw, [:sb, :x, :y])
    push!(puckdfs, raw)
end

# assert that every puck has no missing values and that all spatial barcodes are exactly 14bp and unique
function assert_not_missing(df, col_name, puck)
    @assert all(!ismissing, df[!, col_name]) trouble(
        "puck column $(col_name) in $(puck) contains missing values",
        "Every row must have a barcode and a numeric x and y; a blank or non-numeric cell reads as missing",
        "Columns are positional (spatial barcode, x, y) with no header -- check they are in that order",
        "Inspect the file: `head -3 $(puck)`",
        "Delete the file and re-run to regenerate it from `paths.raw_barcodes_path`",
    )
end
for (puck, puckdf) in zip(pucks, puckdfs)
    logdetail("$(basename(puck)) → $(nrow(puckdf)) spatial barcodes")
    assert_not_missing(puckdf, :sb, puck)
    assert_not_missing(puckdf, :x, puck)
    assert_not_missing(puckdf, :y, puck)
    @assert all(length(s) == 14 for s in puckdf.sb) trouble(
        "some spatial barcodes in $(puck) are not 14bp long",
        "Slide-Tag bead barcodes are 14 bases; a different length means this is not a Slide-Tag puck file",
        "Check the first column really holds barcodes and not, say, an index column",
        "Offending values: $(unique(filter(s -> length(s) != 14, puckdf.sb))[1:min(5, end)])",
        "Delete the file and re-run to regenerate it from `paths.raw_barcodes_path`",
    )
    @assert length(puckdf.sb) == length(Set(puckdf.sb)) trouble(
        "some spatial barcodes in $(puck) are repeated",
        "Each bead barcode must appear once, or its reads cannot be assigned to a single position",
        "A puck CSV concatenated from two sources is the usual cause",
        "Delete the file and re-run to regenerate it from `paths.raw_barcodes_path`",
    )
end

# tag each barcode with its source puck index before concatenating all pucks into one DataFrame
for (i, puckdf) in enumerate(puckdfs)
    puckdf[!, :puck_index] = fill(UInt8(i), nrow(puckdf))
end
puckdf = vcat(puckdfs...)
empty!(puckdfs)  # free memory from per-puck DataFrames

# warn if the same physical barcode sequence appears in more than one puck CSV. This is
# ambiguous: the deduplicated `sb_whitelist` built below (and the SBtoindex matching dict derived
# from it) collapses cross-puck duplicates to a single whitelist entry/index, so a read matching
# that barcode is attributed to only one puck — yet `puckdf` (written out below as puck/sb,
# puck/x, puck/y, puck/puck_index) still keeps one coordinate row per puck occurrence. We do not
# attempt to resolve the ambiguity here (there is no principled way to know which puck a read
# actually came from); we only make the inconsistency visible.
dup_counts = countmap(puckdf.sb)
dup_sbs = [sb for (sb, n) in dup_counts if n > 1]
if !isempty(dup_sbs)
    @warn "$(length(dup_sbs)) spatial barcode(s) appear in more than one puck CSV; reads matching them will be attributed to a single (arbitrary) puck via sb_whitelist/SBtoindex, but the output puck/* coordinate table will still contain one row per puck occurrence for them" example_barcodes=first(dup_sbs, min(5, length(dup_sbs)))
end

# build a sorted, deduplicated whitelist of spatial barcodes
# barcodes with 2+ N bases are excluded as too ambiguous to assign reliably
sb_whitelist = sort(collect(Set{String15}(puckdf.sb)))
num_lowQbeads = count(str -> count(c -> c == 'N', str) >= 2, sb_whitelist)
if num_lowQbeads > 0
    logdetail("dropped $(num_lowQbeads) bead barcode(s) with 2+ N bases")
    sb_whitelist = [str for str in sb_whitelist if count(c -> c == 'N', str) < 2]
end
logts("puck loaded — $(length(sb_whitelist)) unique spatial barcodes")

# UInt32 indices are used throughout; assert the whitelist fits within that range
@assert length(sb_whitelist) < 2^32-1 trouble(
    "$(length(sb_whitelist)) unique spatial barcodes exceeds what a UInt32 index can address",
    "This is a hard limit in spatial_count.jl (sb_i would have to become UInt64), not a problem with your configuration",
    "Report at https://github.com/kabanovskyd/slidr/issues, quoting the barcode count above",
)


##### helper methods ###########################################################

# return the set of all strings within hamming distance `hd` of `str`
function listHDneighbors(str, hd, charlist = ['A','C','G','T','N'])::Set{String}
    res = Set{String}()
    for inds in combinations(1:length(str), hd)
        chars = [str[i] for i in inds]
        pools = [setdiff(charlist, [char]) for char in chars]
        prods = product(pools...)
        for prod in prods
            s = str
            for (i, c) in zip(inds, prod)
                s = s[1:i-1]*string(c)*s[i+1:end]
            end
            push!(res,s)
        end
    end
    return(res)
end

# return all (mutated_string, position) pairs exactly 1 hamming distance from `str`;
# position is the 1-indexed base that was changed
function listHD1neighbors(str, charlist = ['A','C','G','T','N'])::Vector{Tuple{String, UInt8}}
    res = Vector{Tuple{String, UInt8}}()
    for i in 1:length(str)
        for char in setdiff(charlist, str[i])
            s = String15(str[1:i-1]*string(char)*str[i+1:end])
            push!(res, (s,i))
        end
    end
    return(res)
end

# return all strings produced by deleting one character from `str`
function listDEL1neighbors(str)::Vector{String}
    return([string(str[1:i-1], str[i+1:end]) for i in 1:length(str)])
end

# recursively expand all N bases in `s` to A/C/G/T, returning all concrete variants
function expand_N(s::String15, charlist = ['A','C','G','T'])::Vector{String15}
    if !occursin('N', s)
        return [s]
    end
    combins = String15[]
    for nucleotide in charlist
        new_str = String15(replace(s, 'N' => nucleotide, count=1))
        append!(combins, expand_N(new_str))
    end
    return combins
end

# build a lookup dict mapping (sb1[1:8], sb2[9:14]) → (whitelist_index, mismatch_position)
# index == 0 means the match is ambiguous (multiple whitelist entries within fuzzy range)
# index  > 0 means a unique match; pos == 0 for exact, pos != 0 for the corrected base position
#
# also builds a second, fallback dict for single-base deletions in the SECOND half of the
# barcode (sb[9:14]). Such a deletion shifts the downstream sb2 slice left by one base, so the
# observed 6-char sb2 read is [5 true bases][1 unrecoverable base sequenced past the now-shorter
# oligo] — it can never match the primary dict, which is keyed on the full, un-shifted 6-char
# sb2. Instead this fallback is keyed on (sb1, 5-char deletion-variant of sb2), and the caller
# (process_fastqs) tries it against the first 5 of the observed 6 sb2 bases only after the
# primary lookup misses — the 6th observed base is discarded rather than matched against
# anything, since it carries no information about which whitelist barcode this is.
function create_SBtoindex(sb_whitelist)
    SBtoindex = Dict{Tuple{String15, String7}, Tuple{UInt32, Int8}}()
    SBtoindex2 = Dict{Tuple{String15, String7}, Tuple{UInt32, Int8}}()

    # populate fuzzy matches first so exact matches can overwrite them below
    for (i, sb) in enumerate(sb_whitelist)
        numN = count(c -> c == 'N', sb)
        if numN == 0
            # add all hamming-1 neighbours of this barcode
            for (sb_f, ind) in listHD1neighbors(sb)
                sbtuple = (sb_f[1:8], sb_f[9:14])
                haskey(SBtoindex, sbtuple) ? SBtoindex[sbtuple] = (0, 0) : SBtoindex[sbtuple] = (i, ind)
            end
            # add all single-deletion neighbours of the first half of the barcode; mismatch
            # positions here are encoded as -1:-8 (see `l` initialization in process_fastqs)
            sb2 = sb[9:14]
            for (ind, sb1) in enumerate(listDEL1neighbors(sb[1:8]))
                sbtuple = (sb1, sb2)
                haskey(SBtoindex, sbtuple) ? SBtoindex[sbtuple] = (0, 0) : SBtoindex[sbtuple] = (i, -ind)
            end
            # add all single-deletion neighbours of the second half into the fallback dict;
            # mismatch positions here are encoded as -9:-14 to keep them distinct from the
            # first-half deletion range (-1:-8) and the HD1 mismatch-position range (1:14)
            for (ind, sb2_f) in enumerate(listDEL1neighbors(sb2))
                sbtuple2 = (sb[1:8], sb2_f)
                haskey(SBtoindex2, sbtuple2) ? SBtoindex2[sbtuple2] = (0, 0) : SBtoindex2[sbtuple2] = (i, -(8+ind))
            end
        elseif numN == 1
            # barcode has one N — add the N-containing form and all its A/C/G/T expansions
            ind = findfirst(isequal('N'), sb)
            sbtuple = (sb[1:8], sb[9:14])
            haskey(SBtoindex, sbtuple) ? SBtoindex[sbtuple] = (0, 0) : SBtoindex[sbtuple] = (i, ind)
            for sb_f in expand_N(sb)
                sbtuple = (sb_f[1:8], sb_f[9:14])
                haskey(SBtoindex, sbtuple) ? SBtoindex[sbtuple] = (0, 0) : SBtoindex[sbtuple] = (i, ind)
            end
        end
    end

    # exact matches overwrite any fuzzy entry with pos=0, resolving any ambiguity for perfect hits
    for (i, sb) in enumerate(sb_whitelist)
        if !occursin('N', sb)
            sbtuple = (sb[1:8], sb[9:14])
            SBtoindex[sbtuple] = (i, 0)
        end
    end

    return(SBtoindex, SBtoindex2)
end

# encode a 12bp UMI as a UInt32 using 2-bit nucleotide encoding (A=0, C=1, T=2, G=3)
# the px vector holds the positional weights (powers of 4)
const px = [convert(UInt32, 4^i) for i in 0:(12-1)]
function UMItoindex(UMI::StringView{SubArray{UInt8, 1, Vector{UInt8}, Tuple{UnitRange{Int64}}, true}})::UInt32
    return(dot(px, (codeunits(UMI).>>1).&3))
end

# decode a UInt32 back to its 12bp UMI string; bases order must match the encoding above
const bases = ['A','C','T','G'] # MUST NOT change this order
function indextoUMI(i::UInt32)::String15
    return(String15(String([bases[(i>>n)&3+1] for n in 0:2:22])))
end


##### read the FASTQs ##########################################################

# processing steps per read (in order):
#   1. discard reads with an N in the UMI or a homopolymer UMI
#   2. locate the UP site (TCTTCAGCGTTCCCGAGA) in R2 at positions 9:26 with fuzzy matching
#      (allows up to 2 mismatches or 1 deletion in the UP site, or 1 deletion in sb1)
#   3. discard reads whose spatial barcode does not match the puck whitelist
#      (allows 1 mismatch anywhere in the 14bp barcode, or 1 deletion in either half — see
#      create_SBtoindex() for how first- vs second-half deletions are matched differently)
#   4. increment the (cb_i, umi_i, sb_i) count in the sparse count matrix
function process_fastqs(R1s, R2s, sb_whitelist)
    reads = 0
    # m — UP site matching category counts
    m = Dict("exact"=>0, "GG"=>0, "none"=>0, "-1X"=>0, "1D-"=>0, "1D-1X"=>0, "-1D"=>0, "-2X"=>0, "umi_N"=>0, "umi_homopolymer"=>0)
    # p — spatial barcode matching category counts. SB2del/SB2delambig are second-half
    # (sb[9:14]) single-deletion matches/ambiguities, resolved via the SBtoindex2 fallback
    p = Dict("exact"=>0, "HD1"=>0, "HD1ambig"=>0, "SB2del"=>0, "SB2delambig"=>0, "none"=>0)
    # l — position of the corrected base for fuzzy-matched spatial barcodes; -1:-8 are
    # first-half deletion positions, -9:-14 are second-half deletion positions, 1:14 are HD1
    # mismatch positions, 0 is exact
    l = Dict(i=>0 for i in collect(-14:14))
    cb_dictionary = Dict{String31, UInt32}()              # maps cell barcode sequence → integer index
    mat = Dict{Tuple{UInt32, UInt32, UInt32}, UInt32}()   # (cb_i, umi_i, sb_i) → read count

    logts("building barcode matching dictionaries")
    _t_dicts = time()

    # build a set of all UMI sequences within hamming distance 2 of any 12bp homopolymer run
    # these are discarded as likely optical/PCR artefacts
    homopolymer_whitelist = Set{String15}()
    for i in 0:2
        for str in [String15(c^12) for c in ["A","C","G","T","N"]]
            union!(homopolymer_whitelist, listHDneighbors(str, i))
        end
    end

    # pre-compute fuzzy UP site sets so each read only requires set membership tests
    UPseq = String31("TCTTCAGCGTTCCCGAGA")
    UPseqHD1 = Set{String31}(listHDneighbors(UPseq, 1))   # 1 mismatch in UP
    UPseqHD2 = Set{String31}(listHDneighbors(UPseq, 2))   # 2 mismatches in UP
    UPseqLD1 = Set{String31}(listDEL1neighbors(UPseq))    # 1 deletion in UP
    GG_whitelist = Set{String31}(reduce(union, [listHDneighbors("G"^18, i) for i in 0:3]))  # G-run artefacts

    SBtoindex, SBtoindex2 = create_SBtoindex(sb_whitelist)

    logts("✓ dictionaries built — $(fmtduration(time() - _t_dicts))") ; GC.gc()

    logts("reading FASTQs")
    _t_reads = time()

    for fastqpair in zip(R1s, R2s)
        it1 = fastqpair[1] |> open |> GzipDecompressorStream |> FASTQ.Reader;
        it2 = fastqpair[2] |> open |> GzipDecompressorStream |> FASTQ.Reader;
        try
            # manual paired iteration (rather than `for record in zip(it1, it2)`) so that, once
            # the loop ends, we can tell whether it ended because BOTH readers were exhausted
            # together or because one ran out before the other — plain zip() silently truncates
            # to the shorter reader, which would mask a corrupted/truncated mate FASTQ.
            # Advancing both iterators up front (before any `continue` below) is required so that
            # a `continue` never skips advancing past the current record.
            next1 = iterate(it1)
            next2 = iterate(it2)
            while next1 !== nothing && next2 !== nothing
                record1 = next1[1]
                record2 = next2[1]
                next1 = iterate(it1, next1[2])
                next2 = iterate(it2, next2[2])

                reads += 1

                # extract cell barcode (bases 1-16) and UMI (bases 17-28) from R1
                cb  = FASTQ.sequence(record1, 1:16)
                umi = FASTQ.sequence(record1, 17:28)
                if occursin('N', umi)
                    m["umi_N"] += 1
                    continue
                elseif in(umi, homopolymer_whitelist)
                    m["umi_homopolymer"] += 1
                    continue
                end

                # extract the first 32 bases of R2 and locate the UP site to split sb1 and sb2
                # read structure: [sb1 8bp][UP 18bp][sb2 6bp] starting at base 1
                r2 = FASTQ.sequence(record2, 1:32)
                if r2[9:26] == UPseq                    # exact UP match
                    sb1=r2[1:8]; sb2=r2[27:32]; m["exact"]+=1
                elseif in(r2[9:26], UPseqHD1)           # 1 mismatch in UP
                    sb1=r2[1:8]; sb2=r2[27:32]; m["-1X"]+=1
                elseif in(r2[9:26], GG_whitelist)        # G-run artefact — discard
                    m["GG"]+=1; continue
                elseif r2[8:25]==UPseq                   # 1 deletion in sb1, exact UP
                    sb1=r2[1:7]; sb2=r2[26:31]; m["1D-"]+=1
                elseif in(r2[8:25], UPseqHD1)            # 1 deletion in sb1, 1 mismatch in UP
                    sb1=r2[1:7]; sb2=r2[26:31]; m["1D-1X"]+=1
                elseif in(r2[9:25], UPseqLD1)            # 1 deletion within UP
                    sb1=r2[1:8]; sb2=r2[26:31]; m["-1D"]+=1
                elseif in(r2[9:26], UPseqHD2)            # 2 mismatches in UP
                    sb1=r2[1:8]; sb2=r2[27:32]; m["-2X"]+=1
                else                                      # no detectable UP sequence — discard
                    m["none"]+=1; continue
                end

                # look up the spatial barcode against the whitelist
                # sb_i > 0: unique match; sb_i == 0: ambiguous; sb_i == -1: no match
                sb_i, ind = get(SBtoindex, (sb1, sb2), (-1, 0))
                if sb_i > 0 && ind == 0
                    p["exact"]+=1 ; l[ind]+=1
                elseif sb_i > 0
                    p["HD1"]+=1 ; l[ind]+=1
                elseif sb_i == 0
                    p["HD1ambig"]+=1
                    continue
                else
                    # no match against the full 6bp sb2 — try the second-half single-deletion
                    # fallback against the first 5 of the observed 6 sb2 bases (the 6th is
                    # unrecoverable junk sequenced past a now-shorter true oligo)
                    sb_i, ind = get(SBtoindex2, (sb1, sb2[1:5]), (-1, 0))
                    if sb_i > 0
                        p["SB2del"]+=1 ; l[ind]+=1
                    elseif sb_i == 0
                        p["SB2delambig"]+=1
                        continue
                    else
                        p["none"]+=1
                        continue
                    end
                end

                # assign a compact integer index to this cell barcode (first-seen order)
                # then increment the (cb, umi, sb) entry in the sparse count matrix
                cb_i = get!(cb_dictionary, cb, length(cb_dictionary) + 1)
                umi_i = UMItoindex(umi)
                key = (cb_i, umi_i, sb_i)
                mat[key] = get(mat, key, 0) + 1
            end

            if (next1 === nothing) != (next2 === nothing)
                error(trouble(
                    "R1/R2 files have mismatched read counts — possible truncated/corrupted FASTQ ($(fastqpair[1]) / $(fastqpair[2]))",
                    "Paired FASTQs must hold the same number of records; one file ended before the other",
                    "Verify both files decompress cleanly: `gzip -t $(fastqpair[1]) && gzip -t $(fastqpair[2])`",
                    "Count the records if they do: `zcat FILE | wc -l` (should be 4 lines per read, and equal for both)",
                    "Re-copy or re-demultiplex the affected library and re-run with --spatial-count --force",
                ))
            end
        finally
            close(it1)
            close(it2)
        end
    end

    logts("✓ FASTQs read — $(fmtduration(time() - _t_reads)), $(reads) read(s)") ; GC.gc()

    logts("assembling the count matrix")
    _t_matrix = time()

    # flatten the sparse count dict into a DataFrame for serialisation. The final row count
    # (length(mat)) is known ahead of time, so pre-allocate typed column vectors and fill by
    # index rather than growing a DataFrame row-by-row via push! (a DataFrames.jl anti-pattern
    # at scale), then build the DataFrame once at the end.
    n_rows = length(mat)
    cb_i_col  = Vector{UInt32}(undef, n_rows)
    umi_i_col = Vector{UInt32}(undef, n_rows)
    sb_i_col  = Vector{UInt32}(undef, n_rows)
    reads_col = Vector{UInt32}(undef, n_rows)
    for (row, key) in enumerate(keys(mat))
        value = pop!(mat, key)
        cb_i_col[row]  = key[1]
        umi_i_col[row] = key[2]
        sb_i_col[row]  = key[3]
        reads_col[row] = value
    end
    df = DataFrame(cb_i = cb_i_col, umi_i = umi_i_col, sb_i = sb_i_col, reads = reads_col)

    # convert the cb_dictionary to an ordered whitelist DataFrame and verify index contiguity
    cb_whitelist = DataFrame(cb = collect(String31, keys(cb_dictionary)), cb_i = collect(UInt32, values(cb_dictionary)))
    sort!(cb_whitelist, :cb_i)
    @assert cb_whitelist.cb_i == 1:size(cb_whitelist, 1) trouble(
        "the cell-barcode whitelist indices are not contiguous from 1 to $(size(cb_whitelist, 1))",
        "This is an internal indexing bug in spatial_count.jl, not a problem with your data",
        "Report at https://github.com/kabanovskyd/slidr/issues, quoting this whole message",
    )

    logts("✓ matrix assembled — $(fmtduration(time() - _t_matrix)), $(n_rows) row(s)") ; GC.gc()

    # return: total reads, UP filter stats, SB matching stats, mismatch position counts, cb list, count matrix
    return(reads, m, p, l, cb_whitelist.cb, df)
end

reads, m, p, l, cb_whitelist, df = process_fastqs(R1s, R2s, sb_whitelist)

# Verify internal consistency of the counting statistics. Unlike the input checks above, a failure
# here is a bug in this script's own read accounting rather than anything wrong with the user's data,
# so they share one message that says so instead of suggesting fixes that would not help.
const BOOKKEEPING_BUG = [
    "This is an internal read-accounting inconsistency in spatial_count.jl, not a problem with your data",
    "Nothing about your inputs or configuration will fix it -- please report it",
    "Report at https://github.com/kabanovskyd/slidr/issues, quoting this whole message and the counts printed above",
]

@assert reads == sum(values(m)) trouble(
    "read accounting does not balance: $(reads) reads processed but the UP-filter categories sum to $(sum(values(m)))",
    BOOKKEEPING_BUG...,
)
@assert m["1D-"] + m["1D-1X"] + m["exact"] + m["-1X"] + m["-1D"] + m["-2X"] == sum(values(p)) trouble(
    "read accounting does not balance: UP-matched reads sum to $(m["1D-"] + m["1D-1X"] + m["exact"] + m["-1X"] + m["-1D"] + m["-2X"]) but the spatial-barcode categories sum to $(sum(values(p)))",
    BOOKKEEPING_BUG...,
)
@assert p["exact"] + p["HD1"] + p["SB2del"] == sum(df.reads) trouble(
    "read accounting does not balance: $(p["exact"] + p["HD1"] + p["SB2del"]) barcode-matched reads but the count matrix holds $(sum(df.reads))",
    BOOKKEEPING_BUG...,
)
@assert p["exact"] == l[0] && sum(values(l))-l[0] == p["HD1"] + p["SB2del"] trouble(
    "read accounting does not balance: the mismatch-position counts disagree with the spatial-barcode match categories",
    BOOKKEEPING_BUG...,
)
l = sort(DataFrame(k = l |> keys |> collect, v = l |> values |> collect), :k)

# estimate a subsampling curve by simulating random downsampling of the count matrix
# at 21 evenly-spaced fractions from 0% to 100% of total reads
downsampling = UInt32[]
table = countmap(df.reads)
for prob in 0:0.05:1
    s = [length(unique(floor.(sample(0:k*v-1, round(Int,k*v*prob), replace=false)/k))) for (k,v) in zip(keys(table),values(table))]
    append!(downsampling, sum(s))
    GC.gc()
end


##### save results #############################################################

logts("writing SBcounts.h5")
_t_save = time()

outpath = joinpath(puckpath, "SBcounts.h5")
h5open(outpath, "w") do file
    # barcode whitelists used for decoding the integer-indexed matrix
    create_group(file, "lists")
    file["lists/cb_list", compress=9] = cb_whitelist  # Vector{String31}
    file["lists/sb_list", compress=9] = sb_whitelist  # Vector{String15}
    file["lists/puck_list"] = basename.(pucks)         # Vector{String}

    # sparse count matrix stored as parallel index arrays (COO format)
    create_group(file, "matrix")
    file["matrix/cb_index", compress=9] = df.cb_i  # Vector{UInt32}
    file["matrix/umi", compress=9] = df.umi_i      # Vector{UInt32}
    file["matrix/sb_index", compress=9] = df.sb_i  # Vector{UInt32}
    file["matrix/reads", compress=9] = df.reads    # Vector{UInt32}

    # puck barcode coordinates for spatial mapping
    create_group(file, "puck")
    file["puck/sb", compress=9] = puckdf.sb                  # Vector{String15}
    file["puck/x", compress=9] = puckdf.x                    # Vector{Float64}
    file["puck/y", compress=9] = puckdf.y                    # Vector{Float64}
    file["puck/puck_index", compress=9] = puckdf.puck_index  # Vector{UInt8}

    # run-level metadata for provenance and QC
    create_group(file, "metadata")
    file["metadata/R1s"] = R1s
    file["metadata/R2s"] = R2s
    file["metadata/switch"] = convert(Int8, switch)
    file["metadata/num_reads"] = reads
    file["metadata/num_lowQbeads"] = num_lowQbeads

    # UP site matching breakdown (exact, fuzzy, GG artefact, etc.)
    create_group(file, "metadata/UP_matching")
    file["metadata/UP_matching/type"] = keys(m) |> collect
    file["metadata/UP_matching/count"] = values(m) |> collect

    # spatial barcode matching breakdown and per-position mismatch counts
    create_group(file, "metadata/SB_matching")
    file["metadata/SB_matching/type"] = keys(p) |> collect
    file["metadata/SB_matching/count"] = values(p) |> collect
    file["metadata/SB_matching/position"] = l.k |> collect
    file["metadata/SB_matching/position_count"] = l.v |> collect

    file["metadata/downsampling"] = downsampling  # Vector{UInt32}
end;

logts("✓ spatial barcode counting complete — $(fmtduration(time() - _t_save)) writing") ; GC.gc()
