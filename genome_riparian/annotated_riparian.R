# =============================================================================
# annotated_riparian()
# Beautiful DEEPSPACE-style synteny ribbons with per-genome annotation tracks
#
# Dependencies: ggplot2, data.table, GenomicRanges (for GRanges input)
#
# Author: based on DEEPSPACE (Lovell et al.) ribbon geometry, extended here
#         for asymmetric annotation tracks.
# =============================================================================


# =============================================================================
# annotated_riparian.R  —  final clean version
# =============================================================================

utils::globalVariables(c(
  "genome", "chr", "xmin", "xmax", "ymin", "ymax", "y_mid", "x_mid",
  "label", "x", "y", "id", "col", "value_norm", "fill_col",
  "qname", "tname", "qstart", "qend", "qlen", "tstart", "tend", "tlen",
  "query", "target", "blkID", "strand", "nmatch", "alen",
  "ref_chr", "ref_chr_q", "ref_chr_t", "braid_col",
  "bx1q", "bx2q", "bx1t", "bx2t", "qy", "ty", "qxmin", "txmin",
  "grp", "t_mid", "chr_len", "chr_len_bp",
  "win_start", "win_end", "track_col", "stacked_ymax", "stacked_ymin",
  "y_base", "rx_min", "rx_max", "rect_ymin", "rect_ymax",
  "q_plot_width", "q_chr_len", "t_plot_width", "t_chr_len"
))


# =============================================================================
# SECTION 1 — INTERNAL HELPERS
# =============================================================================

# -----------------------------------------------------------------------------
# .make_braid_poly
# -----------------------------------------------------------------------------
.make_braid_poly <- function(x1q, x2q, x1t, x2t,
                             yq, yt,
                             offset_q  = 0.05,
                             offset_t  = 0.05,
                             n_bez     = 80,
                             curvature = 0.7,
                             id        = 1L) {
  
  y_top  <- yq - offset_q
  y_bot  <- yt + offset_t
  y_span <- y_top - y_bot        # always positive
  
  # Asymmetric control points for pronounced S-curve
  cp1_y  <- y_top - curvature * y_span
  cp2_y  <- y_bot + curvature * y_span
  
  t      <- seq(0, 1, length.out = n_bez)
  
  # y: full cubic Bézier
  bez_y  <- (1-t)^3 * y_top   +
    3*(1-t)^2*t * cp1_y +
    3*(1-t)*t^2 * cp2_y +
    t^3 * y_bot
  
  # x: linear interpolation (degenerate cubic, P1=P0, P2=P3)
  lx <- x1q + t * (x1t - x1q)
  rx <- rev(x2q + t * (x2t - x2q))
  
  data.frame(
    x  = c(lx,       rx),
    y  = c(bez_y, rev(bez_y)),
    id = id
  )
}


# -----------------------------------------------------------------------------
# .read_syn_paf
# -----------------------------------------------------------------------------
.read_syn_paf <- function(paf_files, genome_ids) {
  
  message(sprintf("annotated_riparian: found %d PAF file(s):", length(paf_files)))
  for (f in paf_files) message("  ", basename(f))
  
  paf_cols <- c("qname","qlen","qstart","qend","strand",
                "tname","tlen","tstart","tend","nmatch","alen","mapq")
  
  read_paf_file <- function(path) {
    first <- readLines(path, n = 1L)
    if (!length(first) || !nzchar(first)) {
      warning("Empty PAF file: ", path); return(NULL)
    }
    ncols     <- length(strsplit(first, "\t")[[1]])
    col_names <- c(paf_cols,
                   if (ncols > 12) paste0("tag_", seq_len(ncols - 12)))
    tryCatch(
      data.table::fread(path, sep = "\t", header = FALSE,
                        col.names = col_names, fill = TRUE),
      error = function(e) {
        warning("Failed to read PAF: ", path, "\n  ", e$message); NULL
      }
    )
  }
  
  parse_pair <- function(f, genome_ids) {
    base <- basename(f)
    base <- sub("\\.syn\\.paf(\\.gz)?$", "", base)
    base <- sub("\\.paf(\\.gz)?$",       "", base)
    seps <- c("__windows__vs__", "__vs__", "_vs_", "\\.", "__", "_")
    for (sep in seps) {
      parts <- strsplit(base, sep)[[1]]
      for (qi in seq_len(length(parts) - 1)) {
        q_c <- paste(parts[seq_len(qi)],             collapse = "_")
        t_c <- paste(parts[seq(qi+1, length(parts))], collapse = "_")
        if (q_c %in% genome_ids && t_c %in% genome_ids) return(c(q_c, t_c))
        if (t_c %in% genome_ids && q_c %in% genome_ids) return(c(t_c, q_c))
      }
    }
    found <- genome_ids[sapply(genome_ids, grepl, x = base, fixed = TRUE)]
    if (length(found) == 2) return(found)
    warning("Cannot parse genome pair from: ", basename(f))
    c(NA_character_, NA_character_)
  }
  
  file_meta <- do.call(rbind, lapply(paf_files, function(f) {
    p <- parse_pair(f, genome_ids)
    data.frame(query = p[1], target = p[2], path = f,
               stringsAsFactors = FALSE)
  }))
  
  ng    <- length(genome_ids)
  pairs <- data.frame(query  = genome_ids[-ng],
                      target = genome_ids[-1],
                      stringsAsFactors = FALSE)
  
  blk_list <- vector("list", nrow(pairs))
  
  for (i in seq_len(nrow(pairs))) {
    q  <- pairs$query[i]
    tg <- pairs$target[i]
    
    row <- file_meta[
      !is.na(file_meta$query) &
        ((file_meta$query == q  & file_meta$target == tg) |
           (file_meta$query == tg & file_meta$target == q )), ]
    
    if (nrow(row) == 0) {
      warning(sprintf("No PAF for pair %s <-> %s", q, tg)); next
    }
    if (nrow(row) > 1) {
      warning("Multiple PAFs for ", q, "<->", tg, "; using first.")
      row <- row[1, ]
    }
    
    paf <- read_paf_file(row$path)
    if (is.null(paf) || nrow(paf) == 0) {
      warning("Empty PAF: ", row$path); next
    }
    
    swapped <- (row$query == tg)
    if (swapped) {
      paf_out <- data.frame(
        qname  = paf$tname, qlen   = paf$tlen,
        qstart = paf$tstart, qend  = paf$tend,
        strand = paf$strand,
        tname  = paf$qname, tlen   = paf$qlen,
        tstart = paf$qstart, tend  = paf$qend,
        nmatch = paf$nmatch, alen  = paf$alen,
        mapq   = paf$mapq,
        stringsAsFactors = FALSE)
    } else {
      paf_out <- as.data.frame(paf[, paf_cols, with = FALSE])
    }
    
    paf_out$query  <- q
    paf_out$target <- tg
    blk_list[[i]]  <- paf_out
    
    message(sprintf("  %d alignments for %s <-> %s%s",
                    nrow(paf_out), q, tg,
                    if (swapped) " (swapped)" else ""))
  }
  
  result <- do.call(rbind, Filter(Negate(is.null), blk_list))
  if (is.null(result) || nrow(result) == 0)
    stop("No PAF records loaded. Check genomeIDs match filenames.")
  result
}


# -----------------------------------------------------------------------------
# .collapse_to_blocks
# -----------------------------------------------------------------------------
.collapse_to_blocks <- function(paf_dt, gap_threshold = gapThreshold) {
  
  dt <- data.table::as.data.table(paf_dt)
  message(sprintf("  .collapse_to_blocks: %d rows", nrow(dt)))
  
  # Normalize: ensure tstart < tend always
  swap_idx <- which(dt$tstart > dt$tend)
  if (length(swap_idx) > 0) {
    tmp                 <- dt$tstart[swap_idx]
    dt$tstart[swap_idx] <- dt$tend[swap_idx]
    dt$tend[swap_idx]   <- tmp
    message(sprintf("  normalized %d records where tstart > tend",
                    length(swap_idx)))
  }
  
  # ── Try native blkID ──────────────────────────────────────────────────────
  use_blkID <- FALSE
  if ("blkID" %in% names(dt)) {
    n_blk    <- data.table::uniqueN(dt$blkID)
    n_pairs  <- data.table::uniqueN(dt[, paste(qname, tname)])
    pct_na   <- mean(is.na(dt$blkID))
    use_blkID <- (pct_na < 0.05) && (n_blk > n_pairs * 2)
    message(sprintf("  blkID: %d unique, %d chr-pairs, %.1f%% NA -> use=%s",
                    n_blk, n_pairs, 100 * pct_na, use_blkID))
  }
  
  if (use_blkID) {
    collapsed <- dt[, .(
      qstart = min(qstart), qend = max(qend), qlen = qlen[1],
      tstart = min(tstart), tend = max(tend), tlen = tlen[1],
      nmatch = sum(nmatch), alen = sum(alen),
      strand = ifelse(sum(strand=="-") > sum(strand=="+"), "-", "+")
    ), by = .(query, target, qname, tname, blkID)]
    message(sprintf("  -> %d blocks (blkID path)", nrow(collapsed)))
    return(collapsed)
  }
  
  # ── Fallback: collinearity-based block assignment ─────────────────────────
  message(sprintf("  gap-split (threshold = %g bp)", gap_threshold))
  
  data.table::setorder(dt, query, target, qname, tname, strand, qstart)
  
  # Build group key vector outside data.table to avoid return()-inside-:= error
  grp_vec  <- paste(dt$query, dt$target, dt$qname, dt$tname, dt$strand,
                    sep = "__")
  grp_keys <- unique(grp_vec)
  
  blkID_vec <- character(nrow(dt))
  
  for (g in grp_keys) {
    
    idx <- which(grp_vec == g)
    n   <- length(idx)
    
    if (n == 1L) {
      blkID_vec[idx] <- paste0(g, "__1")
      next
    }
    
    qs <- dt$qstart[idx]
    qe <- dt$qend[idx]
    ts <- dt$tstart[idx]
    te <- dt$tend[idx]
    st <- dt$strand[idx][1]
    
    # Query gap: space between end of record i and start of record i+1
    q_gap <- c(Inf, pmax(0, qs[-1] - qe[-n]))
    
    if (st == "+") {
      # "+" strand: tstart should increase with qstart
      # forward gap = tstart[i+1] - tend[i]
      # large positive = hole in target
      # large negative = target jumped backward = non-collinear
      t_gap_val <- c(Inf, ts[-1] - te[-n])
      
    } else {
      # "-" strand: as qstart increases, target midpoint should DECREASE
      t_mid     <- (ts + te) / 2.0
      t_gap_val <- c(Inf, -diff(t_mid))
    }
    
    # New block when query or target gap exceeds threshold,
    # or target moves in the wrong direction
    new_blk <- (q_gap     >  gap_threshold) |
      (t_gap_val >  gap_threshold) |
      (t_gap_val < -gap_threshold)
    
    blkID_vec[idx] <- paste0(g, "__", cumsum(new_blk))
  }
  
  dt[, blkID := blkID_vec]
  
  # ── Collapse ──────────────────────────────────────────────────────────────
  collapsed <- dt[, .(
    qstart = min(qstart), qend = max(qend), qlen = qlen[1],
    tstart = min(tstart), tend = max(tend), tlen = tlen[1],
    nmatch = sum(nmatch), alen = sum(alen),
    strand = strand[1]
  ), by = .(query, target, qname, tname, blkID)]
  
  message(sprintf("  -> %d blocks", nrow(collapsed)))
  collapsed
}


# -----------------------------------------------------------------------------
# .layout_chrs
# -----------------------------------------------------------------------------
.layout_chrs <- function(blk_dt, genome_ids, gap_prop = 0.05) {
  
  chr_lens <- unique(rbind(
    as.data.frame(blk_dt[, .(genome = query,  chr = qname, len = qlen)]),
    as.data.frame(blk_dt[, .(genome = target, chr = tname, len = tlen)])
  ))
  chr_lens <- chr_lens[!duplicated(paste(chr_lens$genome, chr_lens$chr)), ]
  
  max_size <- max(tapply(chr_lens$len, chr_lens$genome, sum))
  
  layout_list <- vector("list", length(genome_ids))
  
  for (i in seq_along(genome_ids)) {
    gid  <- genome_ids[i]
    cg   <- chr_lens[chr_lens$genome == gid, ]
    if (nrow(cg) == 0) next
    
    # Sort chromosomes by number extracted from name
    num_part <- suppressWarnings(as.numeric(gsub("\\D", "", cg$chr)))
    cg       <- cg[order(num_part, cg$chr, na.last = TRUE), ]
    n_chr    <- nrow(cg)
    gap_each <- (max_size * gap_prop) / (n_chr + 1)
    xpos     <- (max_size - (sum(cg$len) + gap_each * (n_chr + 1))) / 2 + gap_each
    
    rows <- vector("list", n_chr)
    for (j in seq_len(n_chr)) {
      rows[[j]] <- data.frame(
        genome  = gid,
        chr     = cg$chr[j],
        chr_len = cg$len[j],
        xmin    = xpos,
        xmax    = xpos + cg$len[j],
        y       = -(i - 1L),
        stringsAsFactors = FALSE)
      xpos <- xpos + cg$len[j] + gap_each
    }
    layout_list[[i]] <- do.call(rbind, rows)
  }
  
  do.call(rbind, layout_list)
}


# -----------------------------------------------------------------------------
# .blocks_to_plot_coords
# -----------------------------------------------------------------------------
.blocks_to_plot_coords <- function(blk_dt, chr_layout) {
  
  blk <- data.table::as.data.table(data.table::copy(blk_dt))
  lay <- data.table::as.data.table(chr_layout)
  
  # Check if zoomed layout (has reg_start column)
  has_reg <- "reg_start" %in% names(lay)
  
  if (has_reg) {
    q_lay <- lay[, .(query  = genome, qname = chr,
                     qxmin = xmin, qy = y,
                     q_reg_start = reg_start)]
    t_lay <- lay[, .(target = genome, tname = chr,
                     txmin = xmin, ty = y,
                     t_reg_start = reg_start)]
  } else {
    q_lay <- lay[, .(query  = genome, qname = chr, qxmin = xmin, qy = y)]
    t_lay <- lay[, .(target = genome, tname = chr, txmin = xmin, ty = y)]
    # Add dummy reg_start = 0 for non-zoomed case
    q_lay[, q_reg_start := 0]
    t_lay[, t_reg_start := 0]
  }
  
  out <- data.table::merge.data.table(blk, q_lay,
                                      by = c("query","qname"), all.x = TRUE)
  out <- data.table::merge.data.table(out,  t_lay,
                                      by = c("target","tname"), all.x = TRUE)
  
  # Plot x = xmin + (bp_coord - reg_start)
  # For non-zoomed: reg_start = 0, so same as before
  out[, `:=`(
    bx1q = qxmin + (qstart - q_reg_start),
    bx2q = qxmin + (qend   - q_reg_start),
    bx1t = txmin + (tstart - t_reg_start),
    bx2t = txmin + (tend   - t_reg_start)
  )]
  
  # "-" strand: swap target corners
  out[strand == "-", `:=`(
    bx1t = txmin + (tend   - t_reg_start),
    bx2t = txmin + (tstart - t_reg_start)
  )]
  
  n_bad <- sum(is.na(out$bx1q) | is.na(out$bx1t))
  if (n_bad > 0)
    warning(n_bad, " blocks have NA coords — chr name mismatch?")
  
  out
}


# -----------------------------------------------------------------------------
# .assign_braid_colors
# -----------------------------------------------------------------------------
.assign_braid_colors <- function(blk_dt, genome_ids, ref_genome, pal) {
  
  dt      <- data.table::as.data.table(blk_dt)
  ref_idx <- which(genome_ids == ref_genome)
  
  ref_chrs <- unique(c(blk_dt$qname[blk_dt$query  == ref_genome],
                       blk_dt$tname[blk_dt$target == ref_genome]))
  chr_map  <- data.table::data.table(
    genome  = ref_genome,
    chr     = ref_chrs,
    ref_chr = ref_chrs)
  
  propagate <- function(from_g, to_g) {
    pb <- dt[(query == from_g & target == to_g) |
               (query == to_g  & target == from_g)]
    if (nrow(pb) == 0) return(NULL)
    fm <- chr_map[genome == from_g]
    if (nrow(fm) == 0) return(NULL)
    
    if (pb$query[1] == from_g) {
      merged <- data.table::merge.data.table(
        pb[, .(chr = qname, to_chr = tname)],
        fm[, .(chr, ref_chr)], by = "chr", all.x = TRUE)
    } else {
      merged <- data.table::merge.data.table(
        pb[, .(chr = tname, to_chr = qname)],
        fm[, .(chr, ref_chr)], by = "chr", all.x = TRUE)
    }
    nm <- merged[!is.na(ref_chr),
                 .(genome = to_g, chr = to_chr, ref_chr)]
    nm[, .(ref_chr = names(sort(table(ref_chr),
                                decreasing = TRUE))[1]),
       by = .(genome, chr)]
  }
  
  if (ref_idx < length(genome_ids))
    for (k in seq(ref_idx, length(genome_ids) - 1)) {
      nw <- propagate(genome_ids[k], genome_ids[k + 1])
      if (!is.null(nw)) chr_map <- unique(rbind(chr_map, nw))
    }
  if (ref_idx > 1)
    for (k in seq(ref_idx, 2)) {
      nw <- propagate(genome_ids[k], genome_ids[k - 1])
      if (!is.null(nw)) chr_map <- unique(rbind(chr_map, nw))
    }
  
  message(sprintf("  Color map: %d entries", nrow(chr_map)))
  
  q_map <- chr_map[, .(query  = genome, qname = chr, ref_chr_q = ref_chr)]
  t_map <- chr_map[, .(target = genome, tname = chr, ref_chr_t = ref_chr)]
  
  dt <- data.table::merge.data.table(dt, q_map,
                                     by = c("query","qname"),  all.x = TRUE)
  dt <- data.table::merge.data.table(dt, t_map,
                                     by = c("target","tname"), all.x = TRUE)
  
  dt[, ref_chr   := data.table::fcoalesce(ref_chr_q, ref_chr_t)]
  dt[, braid_col := pal[ref_chr]]
  dt[is.na(braid_col), braid_col := pal[1]]
  dt
}

# =============================================================================
# .parse_regions  ── parse region specifications
# =============================================================================
# regions can be:
#   1. A named list: list(g1 = list(chr="Chr01_1", start=10e6, end=20e6),
#                         g2 = list(chr="Chr01_2", start=9e6,  end=19e6))
#   2. A data.frame/data.table with columns: genome, chr, start, end
#   3. NULL = no zoom (plot everything)
#
# Returns a data.table with columns: genome, chr, start, end

.parse_regions <- function(regions) {
  
  if (is.null(regions)) return(NULL)
  
  if (is.data.frame(regions) || data.table::is.data.table(regions)) {
    dt <- data.table::as.data.table(regions)
    need <- c("genome","chr","start","end")
    miss <- setdiff(need, names(dt))
    if (length(miss) > 0)
      stop("regions data.frame missing columns: ", paste(miss, collapse=", "))
    return(dt)
  }
  
  if (is.list(regions)) {
    rows <- lapply(names(regions), function(g) {
      r <- regions[[g]]
      if (is.null(r)) return(NULL)
      # allow either list(chr=..., start=..., end=...) per genome
      # or a data.frame of multiple regions per genome
      if (is.data.frame(r)) {
        r$genome <- g
        return(data.table::as.data.table(r))
      }
      data.table::data.table(
        genome = g,
        chr    = r$chr,
        start  = r$start,
        end    = r$end
      )
    })
    return(data.table::rbindlist(Filter(Negate(is.null), rows)))
  }
  
  stop("regions must be NULL, a data.frame, or a named list.")
}


# =============================================================================
# .clip_blocks_to_regions  ── subset and clip PAF blocks to zoom windows
# =============================================================================
# For each genome pair, keep only blocks that overlap the specified regions,
# and clip block coordinates to the region boundaries.

.clip_blocks_to_regions <- function(blk_dt, regions_dt) {
  
  if (is.null(regions_dt) || nrow(regions_dt) == 0) return(blk_dt)
  
  dt  <- data.table::as.data.table(blk_dt)
  reg <- data.table::as.data.table(regions_dt)
  
  # For each block, check whether query AND/OR target overlap a region.
  # We keep a block if EITHER end overlaps a region.
  # We clip coordinates to the region boundaries.
  
  keep_list <- vector("list", nrow(dt))
  
  for (i in seq_len(nrow(dt))) {
    row <- dt[i]
    
    # Find matching region for query side
    q_reg <- reg[genome == row$query  & chr == row$qname &
                   start  <= row$qend   & end >= row$qstart]
    
    # Find matching region for target side
    t_reg <- reg[genome == row$target & chr == row$tname &
                   start  <= row$tend   & end >= row$tstart]
    
    # Keep block if it overlaps a region on EITHER side
    if (nrow(q_reg) == 0 && nrow(t_reg) == 0) next
    
    row_out <- data.table::copy(row)
    
    # Clip query coordinates if region found
    if (nrow(q_reg) > 0) {
      row_out$qstart <- max(row$qstart, q_reg$start[1])
      row_out$qend   <- min(row$qend,   q_reg$end[1])
    }
    
    # Clip target coordinates if region found
    if (nrow(t_reg) > 0) {
      row_out$tstart <- max(row$tstart, t_reg$start[1])
      row_out$tend   <- min(row$tend,   t_reg$end[1])
    }
    
    # Skip if clipping made block invalid
    if (row_out$qstart >= row_out$qend) next
    if (row_out$tstart >= row_out$tend) next
    
    keep_list[[i]] <- row_out
  }
  
  result <- data.table::rbindlist(Filter(Negate(is.null), keep_list))
  message(sprintf("  .clip_blocks_to_regions: %d -> %d blocks after region filter",
                  nrow(dt), nrow(result)))
  result
}


# =============================================================================
# .layout_chrs_zoomed  ── build layout using REGION coordinates as bounds
# =============================================================================
# Instead of using full chromosome lengths, uses the region start/end
# as the xmin/xmax for each chromosome (or the full chromosome if no
# region is specified for that genome/chr).

.layout_chrs_zoomed <- function(blk_dt, genome_ids, regions_dt,
                                gap_prop = 0.01) {
  
  # Get chromosome lengths from PAF as usual
  chr_lens <- unique(rbind(
    as.data.frame(blk_dt[, .(genome = query,  chr = qname, len = qlen)]),
    as.data.frame(blk_dt[, .(genome = target, chr = tname, len = tlen)])
  ))
  chr_lens <- chr_lens[!duplicated(paste(chr_lens$genome, chr_lens$chr)), ]
  
  # If regions specified, restrict to only the chromosomes in regions
  if (!is.null(regions_dt) && nrow(regions_dt) > 0) {
    reg <- as.data.frame(regions_dt)
    
    # Only keep chr/genome combos that appear in regions
    keep <- paste(chr_lens$genome, chr_lens$chr) %in%
      paste(reg$genome, reg$chr)
    chr_lens <- chr_lens[keep, ]
    
    if (nrow(chr_lens) == 0)
      stop(paste0(
        "No chromosomes remain after region filter.\n",
        "  Check that genome/chr names in regions match the PAF files.\n",
        "  regions genomes: ", paste(unique(reg$genome), collapse=", "), "\n",
        "  regions chrs:    ", paste(unique(reg$chr),    collapse=", "), "\n",
        "  PAF genomes:     ",
        paste(unique(c(blk_dt$query, blk_dt$target)), collapse=", ")
      ))
    
    # Replace chr length with region span for layout scaling
    for (ri in seq_len(nrow(reg))) {
      idx <- which(chr_lens$genome == reg$genome[ri] &
                     chr_lens$chr    == reg$chr[ri])
      if (length(idx) > 0)
        chr_lens$len[idx] <- reg$end[ri] - reg$start[ri]
    }
  }
  
  # Only include genomes that have chromosomes remaining
  genome_ids <- genome_ids[genome_ids %in% chr_lens$genome]
  if (length(genome_ids) == 0)
    stop("No genomes remain after region filter.")
  
  max_size <- max(tapply(chr_lens$len, chr_lens$genome, sum))
  
  layout_list <- vector("list", length(genome_ids))
  
  for (i in seq_along(genome_ids)) {
    gid <- genome_ids[i]
    cg  <- chr_lens[chr_lens$genome == gid, ]
    if (nrow(cg) == 0) next
    
    num_part <- suppressWarnings(as.numeric(gsub("\\D","", cg$chr)))
    cg       <- cg[order(num_part, cg$chr, na.last = TRUE), ]
    n_chr    <- nrow(cg)
    gap_each <- (max_size * gap_prop) / (n_chr + 1)
    xpos     <- (max_size - (sum(cg$len) + gap_each*(n_chr+1))) / 2 + gap_each
    
    rows <- vector("list", n_chr)
    for (j in seq_len(n_chr)) {
      
      # Determine xmin/xmax:
      # If a region is specified for this chr, the plot window IS the region.
      # xmin corresponds to region$start, xmax to region$end.
      # Importantly, we store reg_start so that block coordinates can be
      # offset correctly (bp offset = block_coord - reg_start).
      reg_start <- 0L
      reg_end   <- cg$len[j]   # default = full chr
      
      if (!is.null(regions_dt) && nrow(regions_dt) > 0) {
        m <- regions_dt[genome == gid & chr == cg$chr[j]]
        if (nrow(m) > 0) {
          reg_start <- m$start[1]
          reg_end   <- m$end[1]
        }
      }
      
      rows[[j]] <- data.frame(
        genome    = gid,
        chr       = cg$chr[j],
        chr_len   = reg_end - reg_start,  # visible span
        xmin      = xpos,
        xmax      = xpos + (reg_end - reg_start),
        y         = -(i - 1L),
        reg_start = reg_start,            # NEW: offset for coord conversion
        reg_end   = reg_end,
        stringsAsFactors = FALSE
      )
      xpos <- xpos + (reg_end - reg_start) + gap_each
    }
    layout_list[[i]] <- do.call(rbind, rows)
  }
  
  do.call(rbind, layout_list)
}



.assign_type_colors <- function(blk_dt, type_colors) {
  dt <- data.table::as.data.table(blk_dt)
  dt[, braid_col := type_colors[block_type]]
  dt[is.na(braid_col), braid_col := type_colors["SYN"]]
  dt
}



# =============================================================================
# SECTION 2 — MAIN FUNCTION
# =============================================================================

#' Annotated riparian synteny plot
#'
#' @param paf DEEPSPACE result or list(synFiles=..., genomeIDs=...).
#' @param genomeIDs Display order (top to bottom). NULL = use paf order.
#' @param refGenome Colour-anchor genome. Default: first genome.
#' @param gapWidth Inter-chromosome gap fraction (0-1).
#' @param yBuffer Gap between chr bar edge and braid/track (y-units).
#' @param chrHeight Height of chromosome bars (y-units).
#' @param braidOffset Additional offset when no track is present.
#' @param curvature Bézier curve strength (0=straight, 1=very curved).
#' @param chrDensity Optional density data (GRanges or data.frame).
#' @param chrDensityHeight Height of built-in density track.
#' @param densityValueCol Column to use from chrDensity.
#' @param braidColors Colour palette for braids.
#' @param braidAlpha Braid opacity.
#' @param braidOutline Braid outline colour (NA = none).
#' @param chrFill Chromosome bar fill colour.
#' @param chrOutlineColor Chromosome bar outline colour.
#' @param densityFill Built-in density track fill colour.
#' @param densityOutline Built-in density track outline colour.
#' @param n_bez Number of Bézier interpolation points.
#' @param showChrLabels Draw chromosome labels?
#' @param chrLabelSize Chromosome label font size.
#' @param genomeLabelSize Genome label font size.
#' @param bgColor Plot background colour.
#' @return A ggplot object with attributes chr_layout, blk_dt, genomeIDs,
#'   chrHeight, yBuffer, track_registry.
#' @export
annotated_riparian <- function(
    paf                 = NULL,
    blk_dt              = NULL,
    gapThreshold        = 1e5,
    minBlockSize        = 1e5,
    genomeIDs           = NULL,
    refGenome           = NULL,
    regions             = NULL,
    gapWidth            = 0.01,
    yBuffer             = 0.02,
    chrHeight           = 0.05,
    braidOffset         = 0,
    curvature           = 0.7,
    space_above         = 0,
    space_below         = 0,
    chrDensity          = NULL,
    chrDensityHeight    = 0.3,
    densityValueCol     = NULL,
    colorBy             = "chromosome",
    braidColors         = c("#C4645C","#F5915C","#FFC765","#FCF8D5",
                            "#BEF3F9","#66B8FF","#6666FF","#9C63E1","#F4BDFF"),
    typeColors          = c(SYN   = "#DDDDDD",
                            INV   = "#F1B750",
                            TRANS = "#AECD69",
                            DUP   = "#60BDEF"),
    braidAlpha          = 1,
    braidOutline        = NA,
    chrFill             = "white",
    chrOutlineColor     = "grey50",
    densityFill         = "firebrick",
    densityOutline      = NA,
    n_bez               = 80,
    showChrLabels       = TRUE,
    chrLabelSize        = 2,
    genomeLabelSize     = 7,
    genomeLabelColor    = "white",
    bgColor             = "black") {
  
  for (pkg in c("ggplot2","data.table"))
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Package '", pkg, "' required.")
  
  # ── 1. get blk_dt from PAF or directly ────────────────────────────────────
  if (!is.null(blk_dt)) {
    
    message("annotated_riparian: using pre-built blk_dt ...")
    blk_dt <- data.table::as.data.table(blk_dt)
    
    # When colorBy="type", block_type must exist
    if (colorBy == "type" && !"block_type" %in% names(blk_dt))
      stop(paste0(
        "colorBy='type' requires a 'block_type' column in blk_dt.\n",
        "  Use syri_to_blkdt() to build blk_dt from SyRI output."))
    
    # Minimal required columns
    need <- c("query","target","qname","tname",
              "qstart","qend","qlen","tstart","tend","tlen","strand")
    miss <- setdiff(need, names(blk_dt))
    if (length(miss) > 0)
      stop("blk_dt missing columns: ", paste(miss, collapse = ", "))
    
    # Add blkID if missing
    if (!"blkID" %in% names(blk_dt))
      blk_dt[, blkID := paste(query, target, qname, tname,
                              strand, seq_len(.N), sep = "__")]
    
    # genomeIDs from blk_dt if not provided
    if (is.null(genomeIDs)) {
      genomeIDs <- unique(c(blk_dt$query, blk_dt$target))
      message("  genomeIDs inferred from blk_dt: ",
              paste(genomeIDs, collapse = ", "))
    }
    
    # ref_genome
    ref_genome <- if (!is.null(refGenome)) refGenome else genomeIDs[1]
    
  } else {
    
    # Normal PAF workflow
    if (is.null(paf))
      stop("Either 'paf' or 'blk_dt' must be provided.")
    
    if (is.list(paf) && "mapFilePaths" %in% names(paf)) {
      syn_files  <- paf$mapFilePaths$synFile
      gids_all   <- paf$genomeIDs
      ref_genome <- if (!is.null(refGenome)) refGenome else paf$refGenome
    } else if (is.list(paf) && "synFiles" %in% names(paf)) {
      syn_files  <- paf$synFiles
      gids_all   <- paf$genomeIDs
      ref_genome <- if (!is.null(refGenome)) refGenome else gids_all[1]
    } else {
      stop("'paf' must be DEEPSPACE result or list(synFiles=..., genomeIDs=...).")
    }
    
    if (is.null(genomeIDs))  genomeIDs  <- gids_all
    if (is.null(ref_genome)) ref_genome <- genomeIDs[1]
    
    message("annotated_riparian: reading PAF files ...")
    paf_raw <- .read_syn_paf(syn_files, genomeIDs)
    blk_dt  <- .collapse_to_blocks(paf_raw)
  }
  
  if (!ref_genome %in% genomeIDs)
    stop("refGenome '", ref_genome, "' not in genomeIDs.")
  
  # ── min block size filter  ──────
  if (minBlockSize > 0) {
    n_before <- nrow(blk_dt)
    blk_dt   <- blk_dt[
      (qend - qstart) >= minBlockSize &
        (tend - tstart) >= minBlockSize
    ]
    message(sprintf(
      "  minBlockSize filter (%g bp): %d -> %d blocks",
      minBlockSize, n_before, nrow(blk_dt)))
    if (nrow(blk_dt) == 0)
      stop("No blocks remain after minBlockSize filter.")
  }
  
  # ── 2. region filter ───────────────────────────────────────────────────────
  regions_dt <- .parse_regions(regions)
  
  if (!is.null(regions_dt)) {
    message(sprintf("annotated_riparian: zooming to %d region(s) ...",
                    nrow(regions_dt)))
    print(as.data.frame(regions_dt))
    
    blk_dt <- .clip_blocks_to_regions(blk_dt, regions_dt)
    
    if (nrow(blk_dt) == 0)
      stop("No blocks remain after region filter.")
  }
  
  # ── 3. colour assignment ───────────────────────────────────────────────────
  if (colorBy == "type") {
    message("annotated_riparian: colouring braids by synteny type ...")
    blk_dt <- .assign_type_colors(blk_dt, typeColors)
    
  } else {
    message("annotated_riparian: colouring braids by chromosome ...")
    ref_chrs <- sort(unique(c(
      blk_dt$qname[blk_dt$query  == ref_genome],
      blk_dt$tname[blk_dt$target == ref_genome]
    )))
    pal_fn <- grDevices::colorRampPalette(braidColors)
    pal    <- pal_fn(max(length(ref_chrs), 1))
    names(pal) <- ref_chrs
    
    blk_dt <- .assign_braid_colors(blk_dt, genomeIDs, ref_genome, pal)
  }
  
  # ── 4. layout ──────────────────────────────────────────────────────────────
  message("annotated_riparian: computing layout ...")
  
  if (!is.null(regions_dt)) {
    chr_layout <- .layout_chrs_zoomed(blk_dt, genomeIDs,
                                      regions_dt, gap_prop = gapWidth)
  } else {
    chr_layout <- .layout_chrs(blk_dt, genomeIDs, gap_prop = gapWidth)
  }
  
  blk_coords <- .blocks_to_plot_coords(blk_dt, chr_layout)
  complete   <- blk_coords[!is.na(blk_coords$bx1q) &
                             !is.na(blk_coords$bx1t), ]
  
  message(sprintf("  %d / %d blocks with valid coords",
                  nrow(complete), nrow(blk_coords)))
  
  if (nrow(complete) == 0)
    stop("No blocks have valid plot coordinates.")
  
  # ── 5. braid offsets ───────────────────────────────────────────────────────
  ng       <- length(genomeIDs)
  half_chr <- chrHeight / 2
  
  has_builtin_density <- !is.null(chrDensity)
  if (has_builtin_density && space_below == 0)
    space_below <- chrDensityHeight
  
  offset_bot <- setNames(
    rep(half_chr + yBuffer + space_below + braidOffset, ng), genomeIDs)
  offset_top <- setNames(
    rep(half_chr + yBuffer + space_above + braidOffset, ng), genomeIDs)
  
  message(sprintf("  braid offset_top=%.3f  offset_bot=%.3f",
                  offset_top[1], offset_bot[1]))
  
  # ── 6. safety checks ───────────────────────────────────────────────────────
  missing_q <- setdiff(unique(complete$query),  names(offset_bot))
  missing_t <- setdiff(unique(complete$target), names(offset_top))
  if (length(missing_q) > 0)
    stop("offset_bot missing genomes: ", paste(missing_q, collapse = ", "))
  if (length(missing_t) > 0)
    stop("offset_top missing genomes: ", paste(missing_t, collapse = ", "))
  
  # ── 7. braid polygons ──────────────────────────────────────────────────────
  message(sprintf("  building polygons for %d blocks ...", nrow(complete)))
  braid_polys <- vector("list", nrow(complete))
  
  for (i in seq_len(nrow(complete))) {
    bc <- complete[i, ]
    if (anyNA(c(bc$bx1q, bc$bx2q, bc$bx1t, bc$bx2t,
                bc$qy,   bc$ty))) next
    if (abs(as.numeric(bc$bx2q) - as.numeric(bc$bx1q)) < 1 &&
        abs(as.numeric(bc$bx2t) - as.numeric(bc$bx1t)) < 1) next
    
    poly <- .make_braid_poly(
      x1q      = as.numeric(bc$bx1q),
      x2q      = as.numeric(bc$bx2q),
      x1t      = as.numeric(bc$bx1t),
      x2t      = as.numeric(bc$bx2t),
      yq       = as.numeric(bc$qy),
      yt       = as.numeric(bc$ty),
      offset_q = offset_bot[bc$query],
      offset_t = offset_top[bc$target],
      n_bez    = n_bez,
      curvature = curvature,
      id       = i
    )
    poly$col         <- as.character(bc$braid_col)
    braid_polys[[i]] <- poly
  }
  
  valid_polys <- Filter(Negate(is.null), braid_polys)
  message(sprintf("  %d / %d blocks rendered",
                  length(valid_polys), nrow(complete)))
  
  if (length(valid_polys) == 0) {
    warning("No braid polygons built.")
    braid_df <- data.frame(x = numeric(0), y = numeric(0),
                           id = integer(0), col = character(0),
                           stringsAsFactors = FALSE)
  } else {
    braid_df <- do.call(rbind, valid_polys)
  }
  
  # ── 8. chromosome rectangles ───────────────────────────────────────────────
  chr_rect_df <- data.frame(
    genome = chr_layout$genome,
    chr    = chr_layout$chr,
    xmin   = chr_layout$xmin,
    xmax   = chr_layout$xmax,
    ymin   = chr_layout$y - chrHeight / 2,
    ymax   = chr_layout$y + chrHeight / 2,
    y_mid  = chr_layout$y,
    x_mid  = (chr_layout$xmin + chr_layout$xmax) / 2,
    stringsAsFactors = FALSE
  )
  
  # ── 9. y-axis labels ───────────────────────────────────────────────────────
  genome_info <- do.call(rbind, lapply(genomeIDs, function(g) {
    sub  <- chr_layout[chr_layout$genome == g, ]
    lens <- unique(rbind(
      as.data.frame(blk_dt[blk_dt$query  == g, .(chr = qname, len = qlen)]),
      as.data.frame(blk_dt[blk_dt$target == g, .(chr = tname, len = tlen)])
    ))
    lens <- lens[!duplicated(lens$chr), ]
    data.frame(genome = g, y = sub$y[1], total_bp = sum(lens$len),
               stringsAsFactors = FALSE)
  }))
  
  y_breaks <- genome_info$y
  y_labels <- sprintf("%s\n%.1f Mbp",
                      genome_info$genome,
                      genome_info$total_bp / 1e6)
  
  # ── 10. assemble ggplot ────────────────────────────────────────────────────
  message("annotated_riparian: assembling plot ...")
  
  p <- ggplot2::ggplot()
  
  # Braids
  if (nrow(braid_df) > 0) {
    p <- p +
      ggplot2::geom_polygon(
        data      = braid_df,
        mapping   = ggplot2::aes(x = x, y = y, group = id, fill = col),
        alpha     = braidAlpha,
        colour    = braidOutline,
        linewidth = 0
      ) +
      ggplot2::scale_fill_identity()
  }
  
  # Chromosome bars
  p <- p +
    ggplot2::geom_rect(
      data      = chr_rect_df,
      mapping   = ggplot2::aes(xmin = xmin, xmax = xmax,
                               ymin = ymin, ymax = ymax),
      fill      = chrFill,
      colour    = chrOutlineColor,
      linewidth = 0.25
    )
  
  # Chromosome labels
  if (showChrLabels) {
    p <- p +
      ggplot2::geom_text(
        data    = chr_rect_df,
        mapping = ggplot2::aes(x = x_mid, y = y_mid, label = chr),
        size    = chrLabelSize,
        colour  = "black"
      )
  }
  
  # Built-in legacy density track
  if (has_builtin_density) {
    message("annotated_riparian: adding built-in density track ...")
    p <- .add_builtin_density(
      p            = p,
      chr_density  = chrDensity,
      chr_layout   = chr_layout,
      blk_dt       = blk_dt,
      half_chr     = half_chr,
      y_buffer     = yBuffer,
      track_height = chrDensityHeight,
      value_col    = densityValueCol,
      fill_color   = densityFill,
      outline_col  = densityOutline
    )
  }
  
  # Axes and theme
  p <- p +
    ggplot2::scale_y_continuous(
      breaks = y_breaks,
      labels = y_labels,
      expand = ggplot2::expansion(add = c(0.5, 0.5))
    ) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.005, 0.005))
    ) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = bgColor, colour = NA),
      plot.background  = ggplot2::element_rect(fill = bgColor, colour = NA),
      panel.border     = ggplot2::element_blank(),
      panel.grid       = ggplot2::element_blank(),
      axis.ticks       = ggplot2::element_blank(),
      axis.text.x      = ggplot2::element_blank(),
      axis.title       = ggplot2::element_blank(),
      axis.text.y      = ggplot2::element_text(
        colour = genomeLabelColor,
        size   = genomeLabelSize,
        margin = ggplot2::margin(r = 4)),
      plot.margin = ggplot2::margin(10, 10, 10, 10)
    )
  
  # ── 11. attach metadata ────────────────────────────────────────────────────
  attr(p, "chr_layout")     <- chr_layout
  attr(p, "blk_dt")         <- blk_dt
  attr(p, "genomeIDs")      <- genomeIDs
  attr(p, "chrHeight")      <- chrHeight
  attr(p, "yBuffer")        <- yBuffer
  attr(p, "space_above")    <- space_above
  attr(p, "space_below")    <- space_below
  attr(p, "track_registry") <- list()
  
  p
}


# =============================================================================
# SECTION 3 ── DENSITY TRACK HELPERS
# =============================================================================

# -----------------------------------------------------------------------------
# .get_chr_len_col  ── find the length column in chr_sizes regardless of name
# -----------------------------------------------------------------------------
.get_chr_len_col <- function(chr_sizes) {
  candidates <- c("chr_len", "length", "size", "len", "chr_length",
                  "seqlength", "seq_length", "Length", "Size")
  found <- intersect(candidates, names(chr_sizes))
  if (length(found) == 0)
    stop(paste0(
      "chr_sizes has no recognisable length column.\n",
      "  Expected one of: ", paste(candidates, collapse = ", "), "\n",
      "  Actual columns:  ", paste(names(chr_sizes), collapse = ", ")
    ))
  found[1]
}


# -----------------------------------------------------------------------------
# .prep_density_dt  ── validate, normalise, merge chr lengths
# -----------------------------------------------------------------------------
# Returns a data.table with columns:
#   genome, chr, start, end, chr_len_bp, <stack_cols...>
# Each stack_col is normalised 0-1 WITHIN each genome so that the
# tallest bar across all windows reaches exactly track_height.

.prep_density_dt <- function(density_dt, chr_sizes, stack_cols,
                             renormalise = TRUE) {
  # renormalise = FALSE: use raw values directly as fractions (0-1)
  #               TRUE:  rescale so max value across genome = 1
  # For compositional data (gene+TE+CRM+other = 1), use renormalise=FALSE.
  
  dt <- data.table::as.data.table(data.table::copy(density_dt))
  cs <- data.table::as.data.table(data.table::copy(chr_sizes))
  
  need <- c("genome", "chr", "start", "end")
  miss <- setdiff(need, names(dt))
  if (length(miss) > 0)
    stop("density_dt missing columns: ", paste(miss, collapse = ", "),
         "\n  Available: ", paste(names(dt), collapse = ", "))
  
  missing_cols <- setdiff(stack_cols, names(dt))
  if (length(missing_cols) > 0) {
    warning("stack_cols not in density_dt: ",
            paste(missing_cols, collapse = ", "), " — filled with 0.")
    for (col in missing_cols) dt[, (col) := 0]
  }
  
  len_col <- .get_chr_len_col(cs)
  if (!"genome" %in% names(cs))
    stop("chr_sizes must have a 'genome' column.")
  if (!"chr" %in% names(cs)) {
    alt <- intersect(c("chrom","chromosome","seqname","name","seq"), names(cs))
    if (length(alt) == 0)
      stop("chr_sizes has no 'chr' column.")
    data.table::setnames(cs, alt[1], "chr")
  }
  
  cs_sub <- cs[, .(genome, chr, chr_len_bp = get(len_col))]
  if ("chr_len_bp" %in% names(dt)) dt[, chr_len_bp := NULL]
  dt <- data.table::merge.data.table(dt, cs_sub,
                                     by = c("genome","chr"),
                                     all.x = TRUE, sort = FALSE)
  
  n_na <- sum(is.na(dt$chr_len_bp))
  if (n_na > 0)
    warning(n_na, " density windows have no chr_len_bp match.")
  
  # Replace NA values in stack cols with 0
  for (col in stack_cols) {
    dt[is.na(get(col)), (col) := 0]
  }
  
  if (renormalise) {
    # Check if values look compositional (row sums ≈ 1)
    # If so, warn the user
    row_sums <- rowSums(as.matrix(dt[, stack_cols, with = FALSE]))
    if (mean(abs(row_sums - 1) < 0.01, na.rm = TRUE) > 0.8) {
      message(paste0(
        "  .prep_density_dt: values appear compositional (row sums ≈ 1).\n",
        "  Using raw proportions directly (renormalise=FALSE).\n",
        "  Pass renormalise=TRUE explicitly to override."))
      renormalise <- FALSE
    }
  }
  
  if (renormalise) {
    # Scale so the maximum value in each genome = 1
    for (col in stack_cols) {
      dt[, (col) := {
        vals <- get(col)
        mx   <- max(vals, na.rm = TRUE)
        if (mx == 0) rep(0, .N) else vals / mx
      }, by = genome]
    }
  }
  # If renormalise=FALSE: values are used as-is (must already be 0-1)
  
  dt
}



# -----------------------------------------------------------------------------
# .get_y_base  ── compute where each track starts (in plot y-units)
# -----------------------------------------------------------------------------
# Returns a data.table: genome, chr, y_base
# direction = "above" -> track grows toward +y from chr bar top
# direction = "below" -> track grows toward -y from chr bar bottom

.get_y_base <- function(chr_layout, chr_height, y_buffer,
                        direction, registry, track_gap) {
  # NOTE: space_above/space_below are intentionally NOT used here.
  # Tracks always attach directly to the chromosome bar edge.
  # The space parameters only push braids away from the bar.
  
  lay    <- data.table::as.data.table(chr_layout)
  sign_d <- if (direction == "above") +1 else -1
  half   <- chr_height / 2
  
  # Track starts just outside the chromosome bar edge
  y_base_dt <- lay[, .(
    genome = genome,
    chr    = chr,
    y_base = y + sign_d * (half + y_buffer + track_gap)
  )]
  
  # Accumulate height of all previous tracks on the same side
  for (reg in registry) {
    if (reg$direction == direction) {
      y_base_dt[,
                y_base := y_base + sign_d * (reg$track_height + reg$track_gap)
      ]
    }
  }
  
  y_base_dt
}


# -----------------------------------------------------------------------------
# .compute_stacked_rects  ── density windows -> stacked geom_rect rows
# -----------------------------------------------------------------------------
# stack_cols are drawn innermost-first (first col = closest to chr bar).
# Each col occupies a fraction of track_height proportional to its norm value.
# The fractions are computed PER WINDOW so that they sum correctly.
#
# Returns data.table: genome, chr, xmin, xmax, ymin, ymax, fill_col

.compute_stacked_rects <- function(density_dt,
                                   chr_layout,
                                   stack_cols,
                                   stack_colors,
                                   track_height,
                                   y_base_dt,
                                   direction) {
  
  lay <- data.table::as.data.table(chr_layout)
  dt  <- data.table::as.data.table(data.table::copy(density_dt))
  yb  <- data.table::as.data.table(y_base_dt)
  
  sign_d <- if (direction == "above") +1 else -1
  
  # Merge layout info (xmin, xmax, chr_len) into density table
  lay_sub <- lay[, .(genome, chr, xmin, xmax, lay_chr_len = chr_len)]
  dt <- data.table::merge.data.table(
    dt, lay_sub, by = c("genome", "chr"), all.x = TRUE)
  
  # Merge y_base
  dt <- data.table::merge.data.table(
    dt, yb, by = c("genome", "chr"), all.x = TRUE)
  
  # Drop windows with no layout match
  n_nomatch <- sum(is.na(dt$xmin) | is.na(dt$y_base))
  if (n_nomatch > 0) {
    warning(sprintf("%d density windows have no matching chr in layout — skipped.",
                    n_nomatch))
    dt <- dt[!is.na(xmin) & !is.na(y_base)]
  }
  if (nrow(dt) == 0) {
    warning("No density windows remain after merging with layout.")
    return(data.table::data.table(
      genome=character(), chr=character(),
      xmin=numeric(), xmax=numeric(),
      ymin=numeric(), ymax=numeric(),
      fill_col=character()))
  }
  
  # x coordinates: map bp window positions onto plot x-axis
  # chr_len_bp comes from .prep_density_dt merge with chr_sizes
  # lay_chr_len comes from chr_layout (should be identical, but use chr_len_bp
  # for the scale since that's what was used to build the layout)
  dt[, `:=`(
    rx_min = xmin + (start / chr_len_bp) * (xmax - xmin),
    rx_max = xmin + (end   / chr_len_bp) * (xmax - xmin)
  )]
  
  # Build one data.table per stack column, then rbind
  rect_list <- vector("list", length(stack_cols))
  
  for (k in seq_along(stack_cols)) {
    col     <- stack_cols[k]
    col_col <- if (!is.null(stack_colors) && col %in% names(stack_colors))
      stack_colors[[col]] else "grey70"
    
    norm_vals <- dt[[col]]   # already 0-1 from .prep_density_dt
    
    # Cumulative sum of all previous columns for this window
    # = how far from y_base this column's bottom sits
    if (k == 1L) {
      cum_below <- rep(0, nrow(dt))
    } else {
      prev_cols   <- stack_cols[seq_len(k - 1)]
      cum_below   <- rowSums(as.matrix(dt[, prev_cols, with = FALSE]))
    }
    
    y_inner <- dt$y_base + sign_d * cum_below   * track_height
    y_outer <- dt$y_base + sign_d * (cum_below + norm_vals) * track_height
    
    rect_list[[k]] <- data.table::data.table(
      genome   = dt$genome,
      chr      = dt$chr,
      xmin     = dt$rx_min,
      xmax     = dt$rx_max,
      ymin     = pmin(y_inner, y_outer),
      ymax     = pmax(y_inner, y_outer),
      fill_col = col_col
    )
  }
  
  data.table::rbindlist(rect_list)
}


# -----------------------------------------------------------------------------
# .add_builtin_density  ── legacy single-track helper (unchanged interface)
# -----------------------------------------------------------------------------
.add_builtin_density <- function(p, chr_density, chr_layout, blk_dt,
                                 half_chr, y_buffer, track_height,
                                 value_col, fill_color, outline_col) {
  
  if (inherits(chr_density, "GRanges")) {
    df <- as.data.frame(chr_density)
    if (!"genome" %in% names(df)) {
      if (any(grepl("__", as.character(df$seqnames)))) {
        parts  <- strsplit(as.character(df$seqnames), "__")
        df$genome <- sapply(parts, `[`, 1)
        df$chr    <- sapply(parts, `[`, 2)
      } else {
        stop("GRanges seqnames must be 'genome__chr' or carry a 'genome' mcol.")
      }
    } else {
      df$chr <- as.character(df$seqnames)
    }
    meta <- setdiff(names(df),
                    c("seqnames","start","end","width","strand","genome","chr"))
    if (is.null(value_col)) value_col <- meta[1]
    df$value <- df[[value_col]]
  } else {
    df <- as.data.frame(chr_density)
    if (is.null(value_col))
      value_col <- setdiff(names(df), c("genome","chr","start","end"))[1]
    df$value <- df[[value_col]]
  }
  
  # normalise 0-1 per genome
  df <- do.call(rbind, lapply(split(df, df$genome), function(g) {
    rng <- range(g$value, na.rm = TRUE)
    g$value_norm <- if (diff(rng) == 0) 0 else (g$value - rng[1]) / diff(rng)
    g
  }))
  
  lay <- chr_layout
  df  <- merge(df, lay, by = c("genome","chr"), all.x = TRUE)
  
  # chr_len_bp from blk_dt
  cl <- unique(rbind(
    as.data.frame(blk_dt[, .(genome = query,  chr = qname, chr_len_bp = qlen)]),
    as.data.frame(blk_dt[, .(genome = target, chr = tname, chr_len_bp = tlen)])
  ))
  cl  <- cl[!duplicated(paste(cl$genome, cl$chr)), ]
  df  <- merge(df, cl, by = c("genome","chr"), all.x = TRUE)
  df  <- df[!is.na(df$xmin) & !is.na(df$chr_len_bp), ]
  
  df$rx_min    <- df$xmin + (df$start / df$chr_len_bp) * (df$xmax - df$xmin)
  df$rx_max    <- df$xmin + (df$end   / df$chr_len_bp) * (df$xmax - df$xmin)
  df$rect_ymax <- df$y - half_chr - y_buffer
  df$rect_ymin <- df$rect_ymax - df$value_norm * track_height
  
  p + ggplot2::geom_rect(
    data    = df,
    mapping = ggplot2::aes(xmin = rx_min, xmax = rx_max,
                           ymin = rect_ymin, ymax = rect_ymax),
    fill = fill_color, colour = outline_col, linewidth = 0
  )
}


# =============================================================================
# add_density_track()
# =============================================================================

#' Add a stacked density track to a riparian plot
#'
#' @param p ggplot from annotated_riparian().
#' @param density_dt data.table with columns genome, chr, start, end,
#'   and one column per entry in stack_cols.
#' @param chr_sizes data.frame/table with genome, chr, and a length column.
#' @param chr_layout Optional — extracted from p if NULL.
#' @param stack_cols Character vector of column names to stack.
#' @param track_height Total y-height of this track.
#' @param track_position "above" or "below" the chromosome bar.
#' @param stack_colors Named colour vector (names = stack_cols).
#' @param outline_color Bar outline colour (NA = none).
#' @param outline_size Bar outline width.
#' @param track_gap Extra gap between this track and the previous one.
#' @return Augmented ggplot with updated track_registry attribute.
#' @export
add_density_track <- function(
    p,
    density_dt,
    chr_sizes,
    chr_layout     = NULL,
    stack_cols,
    track_height   = 0.3,
    track_position = "above",
    stack_colors   = NULL,
    outline_color  = NA,
    outline_size   = 0,
    track_gap      = 0,
    renormalise    = NULL) {
  
  if (is.null(chr_layout))  chr_layout  <- attr(p, "chr_layout")
  chr_height  <- attr(p, "chrHeight")
  y_buffer    <- attr(p, "yBuffer")
  space_above <- attr(p, "space_above")
  space_below <- attr(p, "space_below")
  registry    <- attr(p, "track_registry")
  
  if (is.null(chr_layout))
    stop("chr_layout not found. Pass explicitly or use annotated_riparian() output.")
  if (is.null(chr_height))  chr_height  <- 0.05
  if (is.null(y_buffer))    y_buffer    <- 0.02
  if (is.null(space_above)) space_above <- 0
  if (is.null(space_below)) space_below <- 0
  if (is.null(registry))    registry    <- list()
  
  direction <- match.arg(track_position, c("above", "below"))
  
  # Auto-detect renormalise
  if (is.null(renormalise)) {
    dt_check <- data.table::as.data.table(density_dt)
    present  <- intersect(stack_cols, names(dt_check))
    if (length(present) > 1) {
      rs <- rowSums(as.matrix(dt_check[, present, with = FALSE]), na.rm = TRUE)
      renormalise <- !(mean(abs(rs - 1) < 0.05, na.rm = TRUE) > 0.8)
      message(sprintf(
        "  add_density_track: auto-detected renormalise=%s (row sums %s 1)",
        renormalise,
        if (renormalise) "!=" else "~"))
    } else {
      renormalise <- TRUE
    }
  }
  
  message(sprintf("add_density_track: preparing '%s' track (%s) ...",
                  paste(stack_cols, collapse = "+"), direction))
  
  density_dt <- .prep_density_dt(density_dt, chr_sizes, stack_cols,
                                 renormalise = renormalise)
  
  # Tracks always attach to chr bar — no space_above/space_below here
  y_base_dt <- .get_y_base(
    chr_layout = chr_layout,
    chr_height = chr_height,
    y_buffer   = y_buffer,
    direction  = direction,
    registry   = registry,
    track_gap  = track_gap
  )
  
  rects <- .compute_stacked_rects(
    density_dt   = density_dt,
    chr_layout   = chr_layout,
    stack_cols   = stack_cols,
    stack_colors = stack_colors,
    track_height = track_height,
    y_base_dt    = y_base_dt,
    direction    = direction
  )
  
  message(sprintf("  %d rectangles for %d stack columns",
                  nrow(rects), length(stack_cols)))
  
  if (nrow(rects) == 0) {
    warning("add_density_track: no rectangles produced.")
    return(p)
  }
  
  p <- p +
    ggplot2::geom_rect(
      data    = as.data.frame(rects),
      mapping = ggplot2::aes(xmin = xmin, xmax = xmax,
                             ymin = ymin, ymax = ymax,
                             fill = fill_col),
      colour    = outline_color,
      linewidth = outline_size
    ) +
    ggplot2::scale_fill_identity()
  
  # Update registry
  registry[[length(registry) + 1]] <- list(
    direction    = direction,
    track_height = track_height,
    track_gap    = track_gap
  )
  
  attr(p, "chr_layout")     <- chr_layout
  attr(p, "chrHeight")      <- chr_height
  attr(p, "yBuffer")        <- y_buffer
  attr(p, "space_above")    <- space_above
  attr(p, "space_below")    <- space_below
  attr(p, "track_registry") <- registry
  
  p
}


# =============================================================================
# add_signal_track — unchanged interface, just calls add_density_track
# =============================================================================

add_signal_track <- function(
    p,
    signal_dt,
    chr_sizes,
    chr_layout     = NULL,
    signal_col,
    track_height   = 0.15,
    track_position = "below",
    signal_color   = "#FF3366",
    outline_color  = NA,
    outline_size   = 0,
    track_gap      = 0) {
  
  add_density_track(
    p              = p,
    density_dt     = signal_dt,
    chr_sizes      = chr_sizes,
    chr_layout     = chr_layout,
    stack_cols     = signal_col,
    track_height   = track_height,
    track_position = track_position,
    stack_colors   = stats::setNames(signal_color, signal_col),
    outline_color  = outline_color,
    outline_size   = outline_size,
    track_gap      = track_gap,
    renormalise    = TRUE   # signal data is not compositional
  )
}


# =============================================================================
# get_riparian_layout()
# =============================================================================

#' Get chromosome layout from a DEEPSPACE result
#' @inheritParams annotated_riparian
#' @return data.frame with genome, chr, chr_len, xmin, xmax, y
#' @export
get_riparian_layout <- function(paf, genomeIDs = NULL, gapWidth = 0.01) {
  
  if (is.list(paf) && "mapFilePaths" %in% names(paf)) {
    syn_files <- paf$mapFilePaths$synFile
    gids_all  <- paf$genomeIDs
  } else {
    syn_files <- paf$synFiles
    gids_all  <- paf$genomeIDs
  }
  if (is.null(genomeIDs)) genomeIDs <- gids_all
  
  paf_raw <- .read_syn_paf(syn_files, genomeIDs)
  blk_dt  <- .collapse_to_blocks(paf_raw)
  .layout_chrs(blk_dt, genomeIDs, gap_prop = gapWidth)
}
