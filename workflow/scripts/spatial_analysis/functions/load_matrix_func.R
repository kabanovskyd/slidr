
# convenience function for reading HDF5 datasets from a provided internal path
read_data_from_h5 <- function(internal_path){return(h5read(sb_path, internal_path))}

# like read_data_from_h5, but for datasets written as 64-bit integers (Julia Int64) that can exceed
# R's 32-bit integer range on large runs (e.g. total read counts). rhdf5's default 64->32-bit
# conversion silently returns NA on overflow instead of erroring, which corrupts downstream sums.
# Only use this for small metadata scalars/vectors -- it forces double precision, which would
# double memory usage if used on the large matrix/* arrays (uint32, ~1B+ rows).
read_data_from_h5_64 <- function(internal_path){return(h5read(sb_path, internal_path, bit64conversion="double"))}


# some 10X library prep kits complement bases 8-9 of the cell barcode; this remaps to the canonical form
# TC<->AG and GT<->CA swaps correct for the complemented encoding used in certain chemistries
remap_10X_CB <- function(vec) {
  stopifnot(class(vec) == "character")
  basemap = setNames(c("AG","TC","CA","GT"), c("TC","AG","GT","CA"))
  stopifnot(substr(vec,8,9) %in% names(basemap))
  ret = paste0(substr(vec,1,7), basemap[substr(vec,8,9)], substr(vec,10,16))
  stopifnot(len(vec) == len(ret))
  stopifnot(nchar(vec) == nchar(ret))
  return(ret)
}



# enumerate all hamming-distance-1 neighbors of a barcode, including substitutions to N
# N neighbors are included so that barcodes containing sequencing ambiguity codes can be matched
listHD1neighbors <- function(input_string) {
  nucleotides <- c('A','C','G','T','N')
  result <- c()
  for (i in 1:nchar(input_string)) {
    for (nuc in nucleotides[nucleotides != substr(input_string, i, i)]) {
      new_string <- paste0(substr(input_string, 1, i-1), nuc, substr(input_string, i+1, nchar(input_string)))
      result <- c(result, new_string)
    }
  }
  return(result)
}



count_umis <- function(df) {
  # return(df %>% group_by(cb_index, sb_index) %>% summarize(umi=n()) %>% arrange(desc(umi)))
  stopifnot(ncol(df) == 4, names(df) == c("cb_index", "umi_2bit", "sb_index", "reads"))
  # sort so that rows with the same (CB, SB) pair are contiguous
  gdf = df %>% filter(reads > 0) %>% select(cb_index, sb_index) %>% arrange(cb_index, sb_index)
  # guard against an empty post-filter dataframe (e.g. puck mismatch or zero valid reads):
  # bnds would be integer(0), and bnds[[1]] below would throw "subscript out of bounds"
  if (nrow(gdf) == 0) {
    gdf$umi = integer(0)
    return(gdf)
  }
  # detect group boundaries: a boundary occurs wherever the next row belongs to a different (CB, SB) pair
  # bnds is the vector of row positions where each (CB, SB) group ends
  bnds = (gdf$cb_index!=lead(gdf$cb_index) | gdf$sb_index!=lead(gdf$sb_index)) %>% tidyr::replace_na(T) %>% which
  # gdf is already sorted, so each group's rows are identical; take the boundary row as the
  # representative instead of dplyr::distinct(), whose vctrs-based hash table can't be built
  # once nrow(gdf) approaches R's 2^31-1 vector-length limit (hits "Error in vec_unique_loc")
  gdf = gdf[bnds, ]
  # umi count per group = distance between consecutive boundary positions; the first group's count is bnds[[1]]
  gdf$umi = (bnds - lag(bnds)) %>% tidyr::replace_na(bnds[[1]])
  gdf %<>% arrange(desc(umi))
  return(gdf)
}


source(paste0(r_func_path, '/plot_utils_func.R'))


fuzzy_matching <- function(df, cb_list, cb_whitelist) {
  # cb_list: all barcodes observed in the SB library (from HDF5 lists/cb_list)
  # cb_whitelist: called-cell barcodes from Cell Ranger (or CellBender), without lane suffix
  stopifnot(!any(duplicated(cb_list)))
  stopifnot(!any(duplicated(cb_whitelist)))
  stopifnot(names(df) == c("cb_index","umi_2bit","sb_index","reads"))

  # detect whether this library uses complemented bases at positions 8-9 by comparing read counts
  # with and without remapping; apply the correction if remapped barcodes match more reads
  cbs = cb_list[df$cb_index]
  reads_noremap = df$reads[cbs %in% cb_whitelist] %>% sum
  reads_remap = df$reads[cbs %in% remap_10X_CB(cb_whitelist)] %>% sum
  remap = reads_remap > reads_noremap ; rm(cbs)
  if (remap) {
    print("Remapping CB whitelist")
    cb_whitelist %<>% remap_10X_CB
  }
  stopifnot(!duplicated(cb_whitelist))

  # exact_dict[i] = whitelist position for cb_list barcode i, or NA if not in whitelist
  exact_dict = match(cb_list, cb_whitelist)

  # build HD1 fuzzy dict: for each whitelist barcode enumerate all HD1 neighbors, then keep only
  # neighbors that appear in cb_list but are not themselves in the whitelist (those are exact matches)
  # barcodes that are HD1-neighbors of multiple whitelist entries are flagged as ambiguous and excluded
  dict = imap(cb_whitelist, ~data.frame(original=.y, neighbor=listHD1neighbors(.x))) %>% {do.call(rbind,.)}
  dict %<>% filter(neighbor %in% cb_list) # remove barcodes that don't exist
  dict %<>% filter(!neighbor %in% cb_whitelist) # remove exact matches
  HD1ambig = unique(dict$neighbor[duplicated(dict$neighbor)])
  dict %<>% filter(!neighbor %in% HD1ambig) # remove ambiguous matches
  stopifnot(!any(duplicated(dict$neighbor)))
  fuzzy_dict = dict$original[match(cb_list, dict$neighbor)]
  HD1ambig_dict = !is.na(match(cb_list, HD1ambig))

  # annotate each row with its match category; at most one of exact/HD1/HD1ambig can be non-NA
  df %<>% mutate(exact = exact_dict[cb_index], HD1 = fuzzy_dict[cb_index], HD1ambig = HD1ambig_dict[cb_index])
  stopifnot((!is.na(df$exact)) + (!is.na(df$HD1)) + df$HD1ambig <= 1)
  rm(dict, exact_dict, fuzzy_dict, HD1ambig_dict) ; invisible(gc())

  # record how many reads fell into each matching category before converting cb_index values
  cb_matching_type = c("exact", "HD1", "HD1ambig","none")
  cb_matching_count = c(df$reads[!is.na(df$exact)] %>% sum,
                        df$reads[!is.na(df$HD1)] %>% sum,
                        df$reads[df$HD1ambig] %>% sum,
                        df$reads[is.na(df$exact) & is.na(df$HD1) & !df$HD1ambig] %>% sum)
  stopifnot(sum(df$reads) == sum(cb_matching_count))

  # convert cb_index to whitelist positions:
  #   df1: exact matches — cb_index becomes the whitelist index
  #   df2: HD1 matches — cb_index becomes the whitelist index; reads are aggregated after merging
  #   df3: unmatched — cb_index is negated to flag as non-whitelist (filtered downstream)
  df1 = df %>% filter(!is.na(exact)) %>% mutate(cb_index = exact) %>% select(1:4)
  df2 = df %>% filter(!is.na(HD1)) %>% mutate(cb_index = HD1) %>% select(1:4)
  df3 = df %>% filter(is.na(exact) & is.na(HD1)) %>% mutate(cb_index = -cb_index) %>% select(1:4)
  df2 %<>% group_by(cb_index, umi_2bit, sb_index) %>% summarize(reads=sum(reads), .groups="drop")
  df12 <- full_join(df1, df2, by = c("cb_index","umi_2bit","sb_index"))
  df12$reads.x %<>% tidyr::replace_na(0) ; df12$reads.y %<>% tidyr::replace_na(0)
  df12 %<>% mutate(reads = reads.x + reads.y) %>% select(-reads.x, -reads.y)
  stopifnot(colnames(df12) == colnames(df3))
  df = rbind(df12, df3)
  stopifnot(df$cb_index != 0)  # zero would be ambiguous between matched and unmatched
  stopifnot(sum(df$reads) == sum(cb_matching_count))

  # update global metadata via <<- (required because metadata is defined in the parent environment)
  meta <- metadata
  meta$SB_info$remap_10X_CB = remap
  meta$CB_matching = setNames(cb_matching_count, cb_matching_type)
  metadata <<- meta

  return(df)
}




remove_chimeras <- function(df) {
  # a chimeric read has the same (CB, UMI) as another read but maps to a different spatial barcode
  # chimeras arise from barcode collision or index hopping; detected by sorting and checking neighbors
  df %<>% arrange(cb_index, umi_2bit, desc(reads))
  # a row is chimeric if its (CB, UMI) pair also appears in the preceding row (before_same)
  # or if its read count ties with the following row sharing the same (CB, UMI) — ambiguous origin (after_same)
  before_same = tidyr::replace_na(df$cb_index==lag(df$cb_index) & df$umi_2bit==lag(df$umi_2bit), FALSE)
  after_same = tidyr::replace_na(df$cb_index==lead(df$cb_index) & df$umi_2bit==lead(df$umi_2bit) & df$reads==lead(df$reads), FALSE)
  chimeric = before_same | after_same
  meta <- metadata ; meta$SB_filtering %<>% c(reads_chimeric=sum(df[chimeric,]$reads)) ; metadata <<- meta
  return(df[!chimeric,])
}



plot_rankplots <- function(df, read_data_from_h5, out_path) {
  gdf = count_umis(df)
  
  # p1
  cb.data = gdf %>% group_by(cb_index) %>% summarize(umi=sum(umi)) %>% arrange(desc(umi)) %>% {mutate(.,index=seq_len(nrow(.)), filter="all cell barcodes")}
  cb.data2 = cb.data %>% filter(cb_index>0) %>% {mutate(.,index=seq_len(nrow(.)), filter="called cell barcodes only")}
  sb_pct_in_called_cells = round(sum(filter(cb.data,cb_index>0)$umi)/sum(cb.data$umi)*100,2) %>% paste0("%")
  p1 = ggplot(mapping=aes(x=index, y=umi,col=filter))+geom_line(data=cb.data)+geom_line(data=cb.data2) +
    scale_x_log10()+scale_y_log10()+theme_bw()+ggtitle("SB UMI per cell")+ylab("SB UMI counts")+xlab("Cell barcodes") +
    theme(legend.position = c(0.05, 0.05), legend.justification = c("left", "bottom"), legend.background = element_blank(), legend.spacing.y = unit(0.1,"lines"), legend.title=element_blank()) +
    annotate("text", x = Inf, y = Inf, label = g("SB UMI in called cells: {sb_pct_in_called_cells}"), hjust = 1.02, vjust = 1.33)
  rm(cb.data, cb.data2) ; invisible(gc())
  
  # p2
  sb.data = gdf %>% group_by(sb_index) %>% summarize(umi=sum(umi)) %>% arrange(desc(umi)) %>% {mutate(.,index=seq_len(nrow(.)),filter="all cell barcodes")}
  sb.data2 = gdf %>% filter(cb_index > 0) %>% group_by(sb_index) %>% summarize(umi=sum(umi)) %>% arrange(desc(umi)) %>% {mutate(.,index=seq_len(nrow(.)),filter="called cell barcodes only")}
  p2 = ggplot(mapping=aes(x=index,y=umi,col=filter))+geom_line(data=sb.data)+geom_line(data=sb.data2)+
    scale_x_log10()+scale_y_log10()+theme_bw()+ggtitle("SB UMI per bead")+ylab("SB UMI counts")+xlab("Beads")+
    theme(legend.position = c(0.05, 0.05), legend.justification = c("left", "bottom"), legend.background = element_blank(), legend.spacing.y = unit(0.1,"lines"), legend.title=element_blank())
  rm(sb.data, sb.data2) ; invisible(gc())
  
  # p3
  x = seq(0, 1, 0.05) * read_data_from_h5("metadata/num_reads")/1000000
  plot.df = data.frame(x=x, y=read_data_from_h5("metadata/downsampling")/1000000)
  p3 = ggplot(plot.df, aes(x=x,y=y)) + geom_point() + theme_bw() + 
    xlab("Millions of reads") + ylab("Millions of filtered SB UMIs") + ggtitle("SB downsampling curve")
  
  # p4
  d = table(df$reads) %>% as.data.frame %>% setNames(c("reads","count")) %>% mutate(reads = as.numeric(levels(reads))[reads])
  d %<>% filter(reads < 10) %>% bind_rows(data.frame(reads = 10, count = sum(d$count[d$reads >= 10])))
  total_reads = prettyNum(read_data_from_h5('metadata/num_reads'), big.mark=',')
  sequencing_saturation = round((1 - nrow(df) / sum(df$reads)) * 100, 2) %>% paste0("%")
  p4 = ggplot(d, aes(x=reads, y=count/1000/1000)) + geom_col() +
    theme_bw() + xlab("Reads per UMI") + ylab("Millions of filtered SB UMIs") + ggtitle("SB read depth") + 
    annotate("text", x = Inf, y = Inf, label = g("sequencing saturation = {sequencing_saturation}\ntotal reads = {total_reads}"), hjust = 1.02, vjust = 1.33) +
    scale_x_continuous(breaks=min(d$reads):max(d$reads), labels=(min(d$reads):max(d$reads)) %>% {ifelse(.==10, "10+", .)})
  
  plot = plot_grid(p1, p2, p3, p4, ncol=2)
  
  suppressMessages(suppressWarnings(make.pdf(plot, file.path(out_path, "SB.pdf"), 7, 8)))
  
  meta <- metadata
  meta$SB_info$UMI_pct_in_called_cells = sb_pct_in_called_cells
  meta$SB_info$sequencing_saturation = sequencing_saturation
  metadata <<- meta
  return(T)
}


load_puckdf <- function(read_data_from_h5) {
  # puck/sb, puck/x, puck/y, puck/puck_index are written by spatial_count.jl from the puck CSV files
  puckdf = data.frame(sb=read_data_from_h5("puck/sb"),
                      x=read_data_from_h5("puck/x"),
                      y=read_data_from_h5("puck/y"),
                      puck_index=as.integer(read_data_from_h5("puck/puck_index")))
  dups = unique(puckdf$sb[duplicated(puckdf$sb)])  # beads with the same barcode on multiple spots
  Ns = puckdf$sb[grepl("N", puckdf$sb)]            # beads with ambiguous bases in their barcode

  # mc = count of the most common nucleotide; lr = length of the longest homopolymer run
  # both are indicators of low-complexity or synthesis-error barcodes
  most_character_count <- function(vec) {
    stopifnot(typeof(vec) == "character")
    degenA = str_count(vec, "A")
    degenC = str_count(vec, "C")
    degenG = str_count(vec, "G")
    degenT = str_count(vec, "T")
    degenN = str_count(vec, "N")
    return(pmax.int(degenA, degenC, degenG, degenT, degenN))
  }
  puckdf$mc = most_character_count(puckdf$sb)
  longest_run_length <- function(vec) {
    stopifnot(typeof(vec) == "character", nchar(vec) > 0)
    ret = stringr::str_extract_all(vec, "(.)\\1*")
    return(map_int(ret, ~max(nchar(.))))
  }
  puckdf$lr = longest_run_length(puckdf$sb)

  # split by puck so quality filters and coordinate transformations are applied per-puck
  puckdfs = map(sort(unique(puckdf$puck_index)), ~filter(puckdf, puck_index==.) %>% select(-puck_index))

  # record pre-filter bead counts before removal
  meta <- metadata
  meta$puck_info = list(puck_name = as.character(read_data_from_h5("lists/puck_list")),
                        num_beads = map_int(puckdfs, nrow))

  # remove duplicated barcodes, barcodes containing N, and low-complexity barcodes
  # thresholds: mc > 10 (>10/14 bases the same nucleotide) or lr > 7 (run of 7+ identical bases)
  meta$bead_info$num_dup = map_int(puckdfs, ~sum(.$sb %in% dups))
  meta$bead_info$num_N = map_int(puckdfs, ~sum(.$sb %in% Ns))
  meta$bead_info$num_degen = map_int(puckdfs, ~sum(.$mc > 10 | .$lr > 7))
  puckdfs %<>% map(~filter(., !sb %in% dups))
  puckdfs %<>% map(~filter(., !sb %in% Ns))
  puckdfs %<>% map(~filter(., mc <= 10, lr <= 7))

  # convert raw bead coordinates to microns using a bead-density-dependent scaling factor:
  #   <150k beads: 0.73 µm/unit (standard Slide-seq puck)
  #   150k–600k beads: 1.46 µm/unit (double-density puck)
  #   >600k beads: 0.645 µm/unit (Slide-tag high-density puck)
  get_scaling_factor <- function(bn) {
    if (bn < 150000) {
      k = 0.73
    } else if (bn < 600000) {
      k = 0.73 * 2
    } else {
      k = 0.645
    }
    return(k)
  }
  meta$puck_info$scaling_factors = map_dbl(meta$puck_info$num_beads, get_scaling_factor)
  puckdfs %<>% map2(meta$puck_info$scaling_factors, ~transmute(.x, sb=sb, x=x*.y, y=y*.y))

  # for multi-puck layouts, line up pucks along the x-axis with no gap between them
  # starts[i] is the x-offset for puck i so pucks are arranged side-by-side
  maxs = map_dbl(puckdfs, ~max(.$x))
  mins = map_dbl(puckdfs, ~min(.$x))
  starts = lag(cumsum(maxs-mins)) %>% tidyr::replace_na(0)
  puckdfs %<>% map2(starts, ~mutate(.x, x = x-min(x)+.y)) # line up from y-axis across
  puckdfs %<>% map(~mutate(., y = y-min(y))) # line up on x-axis

  puckdf = do.call(rbind, puckdfs)
  stopifnot(!any(duplicated(puckdf$sb)))
  # puck_boundaries[i] and puck_boundaries[i+1] delimit the x-range of puck i after layout
  meta$puck_info$puck_boundaries = c(starts, max(puckdf$x))

  # sb_index links puck barcodes to the integer indices used in the SB count matrix
  sb_list = read_data_from_h5("lists/sb_list")
  stopifnot(!any(duplicated(sb_list)))
  puckdf$sb_index = match(puckdf$sb, sb_list)
  puckdf %<>% arrange(sb_index) %>% select(sb_index, x, y)

  metadata <<- meta
  return(puckdf)
}



beadplot <- function(sb.data) {
  ggplot(sb.data, aes(x=x, y=y, col=umi)) +
    rasterize(geom_point(size=0.1, shape=16), dpi=200) +
    coord_fixed(ratio=1) +
    theme_classic() +
    labs(x="x (\u00B5m)", y="y (\u00B5m)") +
    scale_color_viridis(trans="log", option="B", name="UMI") + 
    ggtitle(g("SB UMI per bead"))
}



plot_beadplot <- function(df, puckdf, out_path) {
  sb.data = df %>% count_umis %>% group_by(sb_index) %>% summarize(umi=sum(umi), .groups="drop") %>%
    inner_join(y=puckdf, by="sb_index") %>% arrange(umi)
  p1 = beadplot(sb.data) + ggtitle(g("SB UMI per bead (total)"))
  sb.data = df %>% filter(cb_index>0) %>% count_umis %>% group_by(sb_index) %>% summarize(umi=sum(umi), .groups="drop") %>%
    inner_join(y=puckdf, by="sb_index") %>% arrange(umi)
  p2 = beadplot(sb.data) + ggtitle(g("SB UMI per bead (called cells only)"))
  plot = plot_grid(p1, p2, ncol=1)
  suppressMessages(suppressWarnings(make.pdf(plot, file.path(out_path, "beadplot.pdf"), 7, 8)))
  return(T)
}



plot_bead_filtering <- function(bead_umi, puckdf, umi_cutoff, out_path, percent_umi_filter=1) {
  # bead_umi: total UMI per bead (sb_index), summed across all called cells, before the
  # top-N% filter is applied; puckdf supplies bead (x,y) coordinates for the beadplot panels
  plot.data = bead_umi %>% inner_join(puckdf, by="sb_index")
  p1 = beadplot(plot.data %>% arrange(umi)) + ggtitle("SB UMI per bead (unfiltered)")
  p2 = beadplot(plot.data %>% filter(umi <= umi_cutoff) %>% arrange(umi)) + ggtitle(g("SB UMI per bead (top {percent_umi_filter}% filtered)"))

  p3 = ggplot(bead_umi, aes(x=log10(umi))) +
    geom_histogram(bins=100) +
    geom_vline(xintercept=log10(umi_cutoff), linetype="dashed", color="red") +
    theme_bw() + xlab("log10(SB UMI per bead)") + ylab("Number of beads") +
    ggtitle("Bead UMI distribution") +
    annotate("text", x=log10(umi_cutoff), y=Inf, hjust=-0.05, vjust=1.5, color="red",
             label=g("{100 - percent_umi_filter}th percentile cutoff = {add.commas(round(umi_cutoff))} UMIs"))

  plot = plot_grid(plot_grid(p1, p2, ncol=2), p3, ncol=1, rel_heights=c(0.55, 0.45))
  suppressMessages(suppressWarnings(make.pdf(plot, file.path(out_path, "beadfilter.pdf"), 7, 8)))
  return(T)
}



plot_metrics <- function(metadata, out_path) {
  
  plot.df = list(
    c("Total Reads", metadata$SB_filtering[["reads_total"]] %>% add.commas),
    c("Final UMIs", metadata$SB_filtering[["UMIs_final"]] %>% add.commas),
    c("R1<->R2", metadata$SB_info$switch_R1R2),
    c("Remap 10X CB", metadata$SB_info$remap_10X_CB),
    c("Low Q beads", metadata$bead_info$num_lowQ),
    c(g("High UMI beads (top {metadata$bead_info$percent_umi_filter}%)"), metadata$bead_info$num_highUMI %>% add.commas)
  ) %>% {do.call(rbind,.)} %>% as.data.frame %>% setNames(c("Metric", "Value"))
  p_sp = plot_grid(gdraw("Library information", 13), plot.tab(plot.df), ncol=1, rel_heights=c(0.1,0.6))
  
  header = c("Metric", metadata$puck_info$puck_name %>% str_remove("^Puck_") %>% str_remove("\\.csv$"))
  plot.df = list(c("Beads", metadata$puck_info$num_beads %>% add.commas),
                 c("Filtered beads", Reduce(`+`,metadata$bead_info[c("num_dup","num_N","num_degen")]) %>% add.commas),
                 c("Scaling factor", metadata$puck_info$scaling_factors),
                 c("Final UMIs", metadata$puck_info$umi_final %>% add.commas)
  ) %>% {do.call(rbind,.)} %>% as.data.frame %>% setNames(header)
  p_puck = plot_grid(gdraw("Puck information", 13),
                     plot.tab(plot.df),
                     gdraw("", 13), #spacer
                     ncol=1, rel_heights=c(0.1,0.5,0.1))
  
  UP_matching = metadata$UP_matching
  SB_filtering = metadata$SB_filtering
  
  plot.df = data.frame(a=c("exact","fuzzy", "none", "GG"),
                       b=c(UP_matching[["exact"]],
                           sum(UP_matching[c("1D-","1D-1X","-1X","-1D","-2X")]),
                           UP_matching[["none"]],
                           UP_matching[["GG"]]
                       ) %>% {./sum(.)*100} %>% round(2) %>% paste0("%")
  ) %>% arrange(desc(b)) %>% unname
  p_up = plot_grid(gdraw("UP matching", 13), plot.tab(plot.df), ncol=1, rel_heights=c(0.1,0.4))
  
  # SB matching has 6 categories: fold SB2del into "fuzzy" (alongside HD1) and SB2delambig into
  # "ambig" (alongside HD1ambig), matching the UP-matching multi-category rollup above — otherwise
  # the percentages are computed over a smaller-than-true denominator and misrepresent match quality
  plot.df = data.frame(a=c("exact","fuzzy","none","ambig"),
                       b=c(metadata$SB_matching[["exact"]],
                           sum(metadata$SB_matching[c("HD1","SB2del")]),
                           metadata$SB_matching[["none"]],
                           sum(metadata$SB_matching[c("HD1ambig","SB2delambig")])
                       ) %>% {./sum(.)*100} %>% round(2) %>% unname %>% paste0("%")) %>% unname
  p_sb = plot_grid(gdraw("SB matching", 13), plot.tab(plot.df), ncol=1, rel_heights=c(0.1,0.4))
  
  plot.df = data.frame(a=c("exact","fuzzy","none","ambig"),b=metadata$CB_matching[c("exact","HD1","none","HD1ambig")] %>% {./sum(.)*100} %>% round(2) %>% unname %>% paste0("%")) %>% unname
  p_cb = plot_grid(gdraw("CB matching", 13), plot.tab(plot.df), ncol=1, rel_heights=c(0.1,0.4))
  
  plot.df = list(
    c("Invalid UMI", SB_filtering[["reads_noumi"]]),
    c("No UP", SB_filtering[["reads_noup"]]),
    c("No SB", SB_filtering[["reads_nosb"]]),
    c("Chimeric", SB_filtering[["reads_chimeric"]]),
    c("Invalid SB", SB_filtering[["reads_lowQsb"]]),
    c("Uncalled CB", SB_filtering[["reads_uncalled"]])
  ) %>% {do.call(rbind,.)} %>% as.data.frame %>% setNames(c("Filter","Percent"))
  plot.df$Percent = as.numeric(plot.df$Percent) / (SB_filtering[["reads_total"]] - lag(cumsum(plot.df$Percent),1,0))
  # plot.df[1:4,2] = 1-plot.df[1:4,2] # convert from fraction removed to fraction retained
  plot.df$Percent = round(plot.df$Percent * 100, 2) %>% paste0("%")
  p_filter = plot_grid(gdraw("SB filtering", 13), plot.tab(plot.df), ncol=1, rel_heights=c(0.1,0.7))
  
  plot.df = metadata$SB_fuzzy_position %>% {data.frame(pos=as.numeric(names(.)),count=as.numeric(unname(.)))} %>% filter(pos>0)
  plot.df$pos[plot.df$pos>8] %<>% add(18)
  #plot.df %<>% mutate(count = count/1000/1000)
  p_loc = ggplot(plot.df,aes(x=pos,y=count))+geom_col()+theme_bw() +
    geom_rect(aes(xmin=9, xmax=26, ymin=-Inf, ymax=Inf), fill="grey") +
    annotate("text", x=17.5, y=max(plot.df$count, na.rm=T)*0.1, label="UP Site", color="black") +
    xlab("Spatial barcode base position") + ylab("Fuzzy matches") + ggtitle("Location of spatial barcode fuzzy match")
  
  p_R = list(c("R1s", metadata$SB_info$R1s %>% basename %>% str_remove("\\.fastq\\.gz$") %>% str_remove("_001")),
             c("R2s", metadata$SB_info$R2s %>% basename %>% str_remove("\\.fastq\\.gz$") %>% str_remove("_001"))) %>%
    {do.call(rbind,.)} %>% as.data.frame %>% setNames(NULL) %>% plot.tab
  
  plot = plot_grid(
    gdraw("Spatial library metadata",16),
    plot_grid(p_sp, p_puck, ncol=2, rel_widths = c(0.38,0.62)),
    plot_grid(p_up, p_sb, p_cb, ncol=3),
    plot_grid(p_filter, p_loc, ncol=2, rel_widths = c(0.35,0.65)),
    p_R,
    #gdraw(""), #spacer
    ncol=1,
    rel_heights = c(0.04,0.15,0.11,0.17,0.1)
  )
  
  suppressMessages(suppressWarnings({make.pdf(plot, file.path(out_path, "SBmetrics.pdf"), 7, 8)}))

  return(T)
}