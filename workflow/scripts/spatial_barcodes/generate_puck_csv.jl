# written by Matthew Shebet at the Broad Institute

# activate the project environment so all required packages are available
using Pkg
project_path = ENV["JULIA_PROJECT_PATH"]
Pkg.activate(project_path)
Pkg.precompile()

using DataFrames
using CSV
using Dates


# Progress-logging vocabulary, matching spatial_count.jl and the Python side (helpers.log_ts /
# log_detail / log_section) so every log a run produces reads the same way.
const TS_WIDTH = 10

logts(msg) = (println(Dates.format(now(), "HH:MM:SS") * "  " * msg); flush(stdout))
logdetail(msg) = (println(" "^TS_WIDTH * msg); flush(stdout))

function logsection(title, width=74)
    prefix = "── $title "
    println()
    println(prefix * "─"^max(0, width - length(prefix)))
    flush(stdout)
end

function fmtduration(seconds)
    total = round(Int, seconds)
    hours, rem = divrem(total, 3600)
    minutes, secs = divrem(rem, 60)
    hours > 0 && return "$(hours)h$(lpad(minutes, 2, '0'))m"
    minutes > 0 && return "$(minutes)m$(lpad(secs, 2, '0'))s"
    return "$(secs)s"
end

# Report that one puck could not be processed, with actionable troubleshooting bullets, in the same
# "<message> / Troubleshooting: / • ..." layout the Python side of the pipeline uses.
#
# Every failure in this script is per-puck and non-fatal: the puck is skipped and the batch carries
# on. That makes these messages the *only* warning the user gets -- the run does not stop here, it
# stops later in the spatial-count stage with a "puck file not found" error that says nothing about
# why this puck was skipped. So each skip has to explain itself and what to do about it.
function skip_puck(dir, reason, hints...)
    logts("[WARNING]: skipping puck $(dir) — $(reason)")
    if !isempty(hints)
        logdetail("Troubleshooting:")
        for hint in hints
            logdetail(" • $(hint)")
        end
    end
    logdetail(" • This puck is skipped, but the spatial-count stage will fail later with \"puck file not found\" until it is fixed")
    flush(stdout)
end

# LIB_PUCK_IN  — directory containing raw puck subdirectories (BeadBarcodes.txt + BeadLocations.txt)
# LIB_PUCK_PATH — directory where processed puck CSVs are written
inpath = ENV["LIB_PUCK_IN"]
outpath = ENV["LIB_PUCK_PATH"]

logsection("puck generation")

# collect all subdirectory names in the raw puck input directory
indirs = readdir(inpath)

# derive already-processed puck names by stripping the .csv extension from existing outputs
outdirs = [replace(file, ".csv" => "") for file in readdir(outpath)]

# keep only directories that contain both required barcode files
indirs = filter(indir -> begin
    try
        isfile(joinpath(inpath, indir, "BeadBarcodes.txt")) &&
        isfile(joinpath(inpath, indir, "BeadLocations.txt"))
    catch e
        if e isa Base.IOError
            # skip directories that cannot be read due to permission errors
            skip_puck(indir, "cannot be read ($(e.msg))",
                "Check the read permissions on $(joinpath(inpath, indir))",
                "Check `paths.raw_barcodes_path` points at a directory you can read from this machine")
            false
        else
            rethrow()
        end
    end
end, indirs)

logts("scanning $(inpath)")
logdetail("$(length(indirs)) puck library/libraries present")

# exclude directories whose CSVs have already been generated
indirs = filter(x -> !(x in outdirs), indirs)
logdetail("$(length(indirs)) not yet built")

if length(indirs) == 0
    logts("nothing to build — every puck already has a CSV")
    # A no-op here is normal when every puck has already been built, but it is also what the user
    # sees when the puck they actually need is absent from raw_barcodes_path -- in which case the
    # spatial-count stage fails a moment later with "puck file not found" and nothing explains why.
    # Say up front which case this was.
    logdetail("• Every puck directory under $(inpath) already has a CSV in $(outpath)")
    logdetail("• If the puck you need is missing, check the `Puck ID` metadata column matches a subdirectory name there")
    logdetail("• That subdirectory must contain both BeadBarcodes.txt and BeadLocations.txt")
    exit(0)
end

_t_start = time()
for dir in indirs
    # force a full GC pass before each puck to keep peak memory bounded
    GC.gc(true)

    # Declared before the try, not inside it: a Julia `try` block is its own scope, so a variable
    # first assigned within it is invisible to the matching `catch`. The length-mismatch handler below
    # reports these three counts, and without this the handler itself died with an UndefVarError --
    # turning an intended per-puck skip into a crash that abandoned every remaining puck.
    sbs = String[]
    x = String[]
    y = String[]

    try
        # read barcodes — each row in the file is one barcode string, split across columns by
        # base. Validate the column count is consistent before joining each row into a barcode
        # string: a malformed row (wrong number of fields) would otherwise silently join into a
        # barcode of the wrong length with no error.
        empty!(sbs)
        expected_ncols = -1
        for row in CSV.File(joinpath(inpath, dir, "BeadBarcodes.txt"), header=false, delim=',')
            ncols = length(row)
            if expected_ncols == -1
                expected_ncols = ncols
            end
            if ncols != expected_ncols
                error("BeadBarcodes.txt in $dir has a row with $(ncols) fields, expected $(expected_ncols) (inconsistent barcode length) — refusing to join a possibly-malformed row")
            end
            push!(sbs, join(row))
        end

        # read coordinates — first line is x values, second line is y values, comma-separated.
        # readline() returns "" at EOF rather than raising, so check eof() explicitly before each
        # read: otherwise a truncated file could silently yield an empty coordinate line that,
        # combined with an equally-empty barcode list, would still pass the length-match assert
        # below.
        x, y = open(joinpath(inpath, dir, "BeadLocations.txt"), "r") do file
            eof(file) && error("BeadLocations.txt in $dir is empty or missing the x-coordinate line")
            xline = readline(file)
            eof(file) && error("BeadLocations.txt in $dir is missing the y-coordinate line")
            yline = readline(file)
            (split(xline, ","), split(yline, ","))
        end

        # ensure all three arrays are the same length before zipping into a DataFrame
        @assert length(x) == length(y) == length(sbs)

        # the length-match assert above only guards the token *count* — validate each token is
        # actually a finite number too, so a corrupted coordinate line is caught here (and this
        # puck skipped in isolation, per this loop's design) rather than surfacing later as a much
        # more confusing batch-wide CSV parse failure when spatial_count.jl reads this file back
        # in with typed (String15, Float64, Float64) columns.
        # NOTE: tryparse returns `nothing` for genuinely non-numeric tokens but SUCCEEDS for the
        # literal tokens "NaN"/"Inf"/"-Inf" (returning the corresponding non-finite float), so we
        # must reject both `nothing` and non-finite values — otherwise a garbled coordinate written
        # as the text "NaN" (the exact str(nan) corruption pattern of item 133) slips through as a
        # NaN/Inf bead position. The `isnothing(v) ||` short-circuit keeps `!isfinite` off `nothing`.
        x_parsed = tryparse.(Float64, String.(x))
        y_parsed = tryparse.(Float64, String.(y))
        bad_coord = v -> isnothing(v) || !isfinite(v)
        if any(bad_coord, x_parsed) || any(bad_coord, y_parsed)
            error("BeadLocations.txt in $dir contains non-numeric or non-finite (NaN/Inf) x or y coordinate value(s)")
        end

        # write the combined barcode + coordinate table as a headerless CSV
        df = DataFrame(sb=sbs, x=x, y=y)
        CSV.write(joinpath(outpath, "$(dir).csv"), df, header=false)

    catch e
        if e isa Base.IOError || e isa Base.SystemError
            # skip pucks with unreadable files rather than crashing the whole run
            skip_puck(dir, "its files could not be read ($(e.msg))",
                "Check the read permissions on $(joinpath(inpath, dir))",
                "Check both BeadBarcodes.txt and BeadLocations.txt are present and not still being copied",
                "On a networked filesystem, a transient mount failure also looks like this -- retrying may be enough")
            continue
        elseif e isa AssertionError
            # barcode and coordinate counts disagree — puck file is malformed
            skip_puck(dir, "BeadBarcodes/BeadLocations length mismatch (barcodes=$(length(sbs)), x=$(length(x)), y=$(length(y)))",
                "There must be exactly one x and one y coordinate per bead barcode",
                "Check the two files came from the same puck and neither is truncated",
                "Check BeadLocations.txt has exactly two comma-separated lines (all x values, then all y values)",
                "Re-download the puck's barcode files from the vendor")
            continue
        elseif e isa ErrorException
            # malformed row (wrong field count) or unexpected EOF while reading BeadBarcodes.txt /
            # BeadLocations.txt — skip this puck rather than crashing the whole batch
            skip_puck(dir, e.msg,
                "This usually means one of the two input files is truncated or was written in an unexpected format",
                "Inspect them: `head -2 $(joinpath(inpath, dir, "BeadLocations.txt"))` and `head -3 $(joinpath(inpath, dir, "BeadBarcodes.txt"))`",
                "Re-download the puck's barcode files from the vendor",
                "Or supply a pre-built $(dir).csv in `paths.puck_path` to skip generation entirely")
            continue
        else
            rethrow()
        end
    end
end

logts("✓ puck generation complete — $(fmtduration(time() - _t_start)) for $(length(indirs)) puck(s)")
exit(0)
