# =============================================================================
# read_syri.R
# =============================================================================

read_chrlen <- function(files) {
  # files: named vector, names = genome IDs
  dt_list <- lapply(names(files), function(g) {
    dt <- data.table::fread(files[[g]], header = FALSE,
                            col.names = c("chr","size"))
    dt[, genome := g]
    dt
  })
  data.table::rbindlist(dt_list)[, .(genome, chr, size)]
}


read_syri_output <- function(path, query_genome, target_genome) {
  # SyRI output: 12 tab-separated columns, "-" = missing
  col_names <- c("ref_chr","ref_start","ref_end",
                 "ref_seq","query_seq",
                 "query_chr","query_start","query_end",
                 "unique_id","parent_id","type","extra")
  
  dt <- data.table::fread(path, sep = "\t", header = FALSE,
                          col.names = col_names,
                          fill = TRUE, na.strings = "-")
  
  # Remove comment lines
  dt <- dt[!grepl("^#", ref_chr)]
  
  dt[, `:=`(
    query       = query_genome,
    target      = target_genome,
    ref_start   = as.integer(ref_start),
    ref_end     = as.integer(ref_end),
    query_start = as.integer(query_start),
    query_end   = as.integer(query_end)
  )]
  
  message(sprintf("  %s: %d records", basename(path), nrow(dt)))
  print(dt[, .N, by = type][order(-N)])
  
  dt
}


syri_to_blkdt <- function(syri_files,
                          genome_ids,
                          chr_sizes,
                          include_types = c("SYN","INV","TRANS",
                                            "DUP","INVTR","INVDP",
                                            "SYNAL","INVAL","TRANSAL",
                                            "DUPAL","INVTRAL","INVDPAL"),
                          minBlockSize  = 0) {
  
  # ── Read all files ─────────────────────────────────────────────────────────
  syri_list <- lapply(names(syri_files), function(pair) {
    parts <- strsplit(pair, "__")[[1]]
    if (length(parts) != 2)
      stop("syri_files names must be 'query__target', got: ", pair)
    read_syri_output(syri_files[[pair]], parts[1], parts[2])
  })
  
  syri_dt <- data.table::rbindlist(syri_list, fill = TRUE)
  
  # ── Strategy: use AL types for drawing, parent types for classification ────
  #
  # SyRI hierarchy:
  #   Parent (e.g. TRANS):  has ref coords, query coords = NA → NOT drawable
  #   Child  (e.g. TRANSAL): has BOTH ref and query coords  → drawable
  #
  # For SYN/INV: the parent record itself has both coords → use directly
  # For TRANS/DUP/INVTR/INVDP: use the AL child records
  #
  # Rule:
  #   - Records with NA parent_id AND type in (SYN,INV): use directly
  #   - Records with non-NA parent_id AND type in AL types: use for drawing
  #     but inherit block_type from parent type
  
  all_types <- c("SYN","INV","TRANS","DUP","INVTR","INVDP",
                 "SYNAL","INVAL","TRANSAL","DUPAL","INVTRAL","INVDPAL",
                 "HDR","TDM","CPL","CPG")
  
  dt <- syri_dt[type %in% all_types]
  
  # ── Map every type to a visual block_type ─────────────────────────────────
  dt[, block_type := data.table::fcase(
    type %in% c("SYN","SYNAL","HDR"),           "SYN",
    type %in% c("INV","INVAL"),                  "INV",
    type %in% c("TRANS","TRANSAL"),              "TRANS",
    type %in% c("DUP","DUPAL","TDM"),            "DUP",
    type %in% c("INVTR","INVTRAL"),              "INV",    # inverted translocation → INV color
    type %in% c("INVDP","INVDPAL"),              "DUP",    # inverted duplication → DUP color
    type %in% c("CPL","CPG"),                    "TRANS",  # copy-paste → TRANS color
    default = "SYN"
  )]
  
  # ── Select drawable records ────────────────────────────────────────────────
  # A record is drawable if it has BOTH ref and query coordinates
  dt_draw <- dt[
    !is.na(ref_chr)   & !is.na(query_chr) &
      !is.na(ref_start) & !is.na(query_start) &
      ref_start > 0     & query_start > 0
  ]
  
  # For top-level types that carry their own coords (SYN, INV etc.):
  # keep only top-level (parent_id == NA) to avoid double-drawing with AL
  # For AL types: keep all (they are children, parent has no query coords)
  al_types     <- c("SYNAL","INVAL","TRANSAL","DUPAL","INVTRAL","INVDPAL")
  parent_types <- c("SYN","INV","TRANS","DUP","INVTR","INVDP","HDR","TDM","CPL","CPG")
  
  dt_parents <- dt_draw[type %in% parent_types & is.na(parent_id)]
  dt_al      <- dt_draw[type %in% al_types]
  
  # For parent types that have AL children, the parent's query coords are NA
  # so they were already excluded by the !is.na(query_chr) filter above.
  # Just combine:
  dt_final <- data.table::rbindlist(list(dt_parents, dt_al), fill = TRUE)
  
  message(sprintf("\nDrawable blocks: %d", nrow(dt_final)))
  print(dt_final[, .N, by = .(type, block_type)][order(block_type, -N)])
  
  if (nrow(dt_final) == 0)
    stop("No drawable records found.")
  
  # ── Build blk_dt ──────────────────────────────────────────────────────────
  blk_dt <- dt_final[, .(
    query      = query,
    target     = target,
    qname      = ref_chr,
    qstart     = as.integer(ref_start),
    qend       = as.integer(ref_end),
    tname      = query_chr,
    tstart     = as.integer(query_start),
    tend       = as.integer(query_end),
    blkID      = unique_id,
    strand     = data.table::fcase(
      type %in% c("INV","INVAL","INVTR","INVTRAL",
                  "INVDP","INVDPAL"),  "-",
      default = "+"
    ),
    block_type = block_type,
    nmatch     = 0L,
    alen       = as.integer(abs(ref_end - ref_start))
  )]
  
  # Normalize tstart < tend
  swap_idx <- which(blk_dt$tstart > blk_dt$tend)
  if (length(swap_idx) > 0) {
    tmp                      <- blk_dt$tstart[swap_idx]
    blk_dt$tstart[swap_idx] <- blk_dt$tend[swap_idx]
    blk_dt$tend[swap_idx]   <- tmp
    message(sprintf("  Swapped %d records where tstart > tend", length(swap_idx)))
  }
  
  # ── Merge chromosome lengths ───────────────────────────────────────────────
  cs      <- data.table::as.data.table(chr_sizes)
  len_col <- .get_chr_len_col(cs)
  cs_sub  <- cs[, .(genome, chr, size = get(len_col))]
  
  q_len <- cs_sub[, .(query  = genome, qname = chr, qlen = size)]
  t_len <- cs_sub[, .(target = genome, tname = chr, tlen = size)]
  
  blk_dt <- data.table::merge.data.table(blk_dt, q_len,
                                         by = c("query","qname"), all.x = TRUE)
  blk_dt <- data.table::merge.data.table(blk_dt, t_len,
                                         by = c("target","tname"), all.x = TRUE)
  
  blk_dt[is.na(qlen), qlen := max(qend, na.rm=TRUE), by = .(query,  qname)]
  blk_dt[is.na(tlen), tlen := max(tend, na.rm=TRUE), by = .(target, tname)]
  
  message(sprintf("\nsyri_to_blkdt: %d final blocks", nrow(blk_dt)))
  print(blk_dt[, .N, by = block_type][order(-N)])
  
  blk_dt
  
  # ── Filter by minimum block size ──────────────────────────────────────────
  if (minBlockSize > 0) {
    n_before <- nrow(blk_dt)
    blk_dt <- blk_dt[
      (qend - qstart) >= minBlockSize &
        (tend - tstart) >= minBlockSize
    ]
    message(sprintf(
      "  minBlockSize filter (%g bp): %d -> %d blocks removed %d",
      minBlockSize, n_before, nrow(blk_dt), n_before - nrow(blk_dt)))
    print(blk_dt[, .N, by = block_type][order(-N)])
  }
  
  blk_dt
}
