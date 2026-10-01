# =============================================================================
# density_compute.R
#
# Stage 1 — raw density/signal computation
#   calc_windowed_density()   BED features  → windowed coverage or count
#   calc_windowed_signal()    bedGraph/wig  → windowed mean signal
#
# Stage 2 — normalization
#   normalize_density()       min-max per chromosome per feature
#
# Stage 3 — convenience wrappers
#   prepare_feature_density() BED  → compute + normalize in one call
#   prepare_signal_density()  bedGraph → compute + normalize in one call
#
# Normalization note:
#   Only min-max is offered:  (x - min) / (max - min)  per chromosome.
#   This is honest, simple, and the axis label "relative density" is
#   defensible. Values are always [0,1] after normalization.
#   Raw proportions (coverage mode, with "other") sum to 1 and are used
#   directly for stacked plots without normalization.
#
# Dependencies: data.table
# =============================================================================

utils::globalVariables(c(
  "win_id", "win_size_actual", "fs", "fe",
  "overlap_bp", "feat_val", "chr_len_bp", "size"
))


# =============================================================================
# .get_chr_len_col()   internal — find length column in chr_sizes
# =============================================================================
.get_chr_len_col <- function(chr_sizes) {
  candidates <- c("size","chr_len","length","len","chr_length",
                  "seqlength","seq_length","Length","Size")
  found <- intersect(candidates, names(chr_sizes))
  if (length(found) == 0)
    stop("chr_sizes has no recognisable length column.\n",
         "  Expected one of: ", paste(candidates, collapse=", "), "\n",
         "  Got: ", paste(names(chr_sizes), collapse=", "))
  found[1]
}


# =============================================================================
# .build_windows()   internal — tile chromosomes into fixed windows
# =============================================================================
.build_windows <- function(chr_sizes, win_size, step_size) {
  
  cs  <- data.table::as.data.table(chr_sizes)
  lc  <- .get_chr_len_col(cs)
  
  if (!"chr" %in% names(cs)) {
    alt <- intersect(c("chrom","chromosome","seqname","name"), names(cs))
    if (length(alt) == 0) stop("chr_sizes has no 'chr' column.")
    data.table::setnames(cs, alt[1], "chr")
  }
  
  win_list <- vector("list", nrow(cs))
  for (i in seq_len(nrow(cs))) {
    sz  <- as.numeric(cs[[lc]][i])
    ws  <- seq(1, sz, by = step_size)
    we  <- pmin(ws + win_size - 1, sz)
    win_list[[i]] <- data.table::data.table(
      genome          = cs$genome[i],
      chr             = cs$chr[i],
      start           = ws,
      end             = we,
      win_size_actual = we - ws + 1,
      chr_len_bp      = sz
    )
  }
  
  wins <- data.table::rbindlist(win_list)
  wins[, win_id := .I]
  wins
}


# =============================================================================
# calc_windowed_density()
# =============================================================================

#' Compute windowed feature densities from a BED data frame
#'
#' Tiles each chromosome into fixed windows and computes per-window
#' coverage (bp proportion) or feature count for each category.
#'
#' @param features   data.table: genome | chr | start | end | category
#' @param chr_sizes  data.table: genome | chr | size  (bp)
#' @param win_size   Window size in bp. Default 1e5.
#' @param step_size  Step in bp. NULL = non-overlapping.
#' @param categories Categories to quantify. Include "other" to add remainder
#'   so that coverage proportions sum to 1 per window (for stack mode).
#' @param method     "coverage": bp overlap / window bp  → [0,1] proportion.
#'                   "count":    number of features per window → raw integer.
#'
#' @return data.table:
#'   genome | chr | start | end | chr_len_bp | <one col per category>
#'   Coverage mode: values in [0,1], row sums = 1 when "other" included.
#'   Count mode: raw integers.
#' @export
calc_windowed_density <- function(
    features,
    chr_sizes,
    win_size   = 1e5,
    step_size  = NULL,
    categories = c("gene", "CRM", "LINE", "other"),
    method     = c("coverage", "count")
) {
  
  method    <- match.arg(method)
  win_size  <- as.numeric(win_size)
  step_size <- as.numeric(if (is.null(step_size)) win_size else step_size)
  
  dt <- data.table::as.data.table(data.table::copy(features))
  dt[, start := as.numeric(start)]
  dt[, end   := as.numeric(end)]
  
  # validate
  need <- c("genome","chr","start","end","category")
  miss <- setdiff(need, names(dt))
  if (length(miss) > 0)
    stop("features missing columns: ", paste(miss, collapse=", "))
  
  real_cats <- setdiff(categories, "other")
  do_other  <- "other" %in% categories
  
  not_found <- setdiff(real_cats, unique(dt$category))
  if (length(not_found) > 0)
    message("Note: categories not in features (will be 0): ",
            paste(not_found, collapse=", "))
  
  # build windows
  message(sprintf(
    "calc_windowed_density: window=%g bp, step=%g bp, method=%s",
    win_size, step_size, method))
  
  wins <- .build_windows(chr_sizes, win_size, step_size)
  for (cat in categories) wins[, (cat) := 0]
  
  message(sprintf("  %d windows across %d chr/genome combos",
                  nrow(wins),
                  data.table::uniqueN(wins[, paste(genome, chr)])))
  
  # two-pointer coverage sweep
  .sweep_cov <- function(fs, fe, ws, we) {
    nw <- length(ws); nf <- length(fs)
    out <- numeric(nw)
    if (nf == 0L) return(out)
    fi <- 1L
    for (w in seq_len(nw)) {
      while (fi <= nf && fe[fi] < ws[w]) fi <- fi + 1L
      if (fi > nf) break
      j <- fi; tot <- 0
      while (j <= nf && fs[j] <= we[w]) {
        tot <- tot + pmin(fe[j], we[w]) - pmax(fs[j], ws[w]) + 1
        j   <- j + 1L
      }
      out[w] <- tot
    }
    out
  }
  
  .sweep_cnt <- function(fs, fe, ws, we) {
    mids <- (fs + fe) / 2
    vapply(seq_along(ws),
           function(w) sum(mids >= ws[w] & mids <= we[w]),
           numeric(1))
  }
  
  chr_keys <- unique(wins[, .(genome, chr)])
  
  for (cat in real_cats) {
    message(sprintf("  '%s' ...", cat))
    f_cat <- dt[category == cat,
                .(genome, chr, fs=as.numeric(start), fe=as.numeric(end))]
    if (nrow(f_cat) == 0L) next
    
    for (ki in seq_len(nrow(chr_keys))) {
      g  <- chr_keys$genome[ki]; ch <- chr_keys$chr[ki]
      widx <- wins[genome==g & chr==ch, win_id]
      if (length(widx) == 0L) next
      fw <- wins[win_id %in% widx]
      ff <- f_cat[genome==g & chr==ch]
      if (nrow(ff) == 0L) next
      ff <- ff[order(fs)]
      
      vals <- if (method=="coverage")
        .sweep_cov(ff$fs, ff$fe, fw$start, fw$end)
      else
        .sweep_cnt(ff$fs, ff$fe, fw$start, fw$end)
      
      wins[win_id %in% widx, (cat) := vals]
    }
  }
  
  # other = remainder
  if (do_other) {
    message("  'other' = remainder ...")
    rs <- rowSums(wins[, real_cats, with=FALSE])
    wins[, other := pmax(0, win_size_actual - rs)]
  }
  
  # convert to proportions (coverage only)
  if (method == "coverage") {
    for (cat in categories)
      wins[, (cat) := get(cat) / win_size_actual]
  }
  
  wins[, win_size_actual := NULL]
  wins[, win_id          := NULL]
  
  message(sprintf("calc_windowed_density: done. %d windows.", nrow(wins)))
  wins
}


# =============================================================================
# calc_windowed_signal()
# =============================================================================

#' Compute windowed mean signal from a bedGraph or pre-binned table
#'
#' Aggregates continuous signal values (e.g. CenH3 ChIP-seq bedGraph,
#' read depth, methylation) into fixed-size windows by computing the
#' coverage-weighted mean per window.
#'
#' @section Input formats:
#' \describe{
#'   \item{bedGraph}{
#'     Standard bedGraph: genome | chr | start | end | score
#'     Each interval carries a single score value.
#'     The weighted mean per window is: sum(score_i * overlap_i) / window_bp.
#'   }
#'   \item{Pre-binned table}{
#'     Any table with genome | chr | start | end | <signal_col>.
#'     If bins are already the same size as win_size, the output is identical
#'     to the input (just with chr_len_bp added).
#'   }
#' }
#'
#' @param signal_df  data.frame/table:
#'   genome | chr | start | end | <signal_col>
#' @param chr_sizes  data.table: genome | chr | size (bp)
#' @param signal_col Name of the score/signal column. Default "score".
#' @param win_size   Window size in bp. Default 1e5.
#' @param step_size  Step size. NULL = non-overlapping.
#' @param agg_fun    Aggregation: "weighted_mean" (default, coverage-weighted),
#'                   "mean" (simple mean of overlapping intervals),
#'                   "max", "sum".
#'
#' @return data.table:
#'   genome | chr | start | end | chr_len_bp | signal
#'   "signal" column contains the aggregated value (NOT yet normalized).
#' @export
calc_windowed_signal <- function(
    signal_df,
    chr_sizes,
    signal_col = "score",
    win_size   = 1e5,
    step_size  = NULL,
    agg_fun    = c("weighted_mean", "mean", "max", "sum")
) {
  
  agg_fun   <- match.arg(agg_fun)
  win_size  <- as.numeric(win_size)
  step_size <- as.numeric(if (is.null(step_size)) win_size else step_size)
  
  dt <- data.table::as.data.table(data.table::copy(signal_df))
  dt[, start := as.numeric(start)]
  dt[, end   := as.numeric(end)]
  
  # flexible column name
  if (!signal_col %in% names(dt)) {
    alt <- setdiff(names(dt), c("genome","chr","start","end",
                                "chrom","chromosome"))
    if (length(alt) == 0)
      stop("Cannot find signal column '", signal_col,
           "'. Available: ", paste(names(dt), collapse=", "))
    message(sprintf("calc_windowed_signal: '%s' not found, using '%s'",
                    signal_col, alt[1]))
    signal_col <- alt[1]
  }
  dt[, signal_val := as.numeric(get(signal_col))]
  
  need <- c("genome","chr","start","end")
  miss <- setdiff(need, names(dt))
  if (length(miss) > 0)
    stop("signal_df missing columns: ", paste(miss, collapse=", "))
  
  message(sprintf(
    "calc_windowed_signal: column='%s', window=%g bp, agg='%s'",
    signal_col, win_size, agg_fun))
  
  wins <- .build_windows(chr_sizes, win_size, step_size)
  wins[, signal := 0]
  
  chr_keys <- unique(wins[, .(genome, chr)])
  
  for (ki in seq_len(nrow(chr_keys))) {
    g  <- chr_keys$genome[ki]; ch <- chr_keys$chr[ki]
    widx <- wins[genome==g & chr==ch, win_id]
    if (length(widx) == 0L) next
    fw <- wins[win_id %in% widx]
    ff <- dt[genome==g & chr==ch]
    if (nrow(ff) == 0L) next
    
    sig_vec <- numeric(nrow(fw))
    
    for (w in seq_len(nrow(fw))) {
      ws <- fw$start[w]; we <- fw$end[w]
      
      # overlapping intervals
      ov <- ff[start < we & end > ws]
      if (nrow(ov) == 0L) next
      
      if (agg_fun == "weighted_mean") {
        # weight by overlap length
        ol  <- pmin(ov$end, we) - pmax(ov$start, ws)
        ol  <- pmax(ol, 0)
        tot <- sum(ol)
        sig_vec[w] <- if (tot > 0) sum(ov$signal_val * ol) / tot else 0
        
      } else if (agg_fun == "mean") {
        sig_vec[w] <- mean(ov$signal_val, na.rm=TRUE)
        
      } else if (agg_fun == "max") {
        sig_vec[w] <- max(ov$signal_val, na.rm=TRUE)
        
      } else {
        sig_vec[w] <- sum(ov$signal_val, na.rm=TRUE)
      }
    }
    
    wins[win_id %in% widx, signal := sig_vec]
  }
  
  # clean up
  wins[, win_size_actual := NULL]
  wins[, win_id          := NULL]
  data.table::setnames(wins, "signal", signal_col)
  
  message(sprintf("calc_windowed_signal: done. %d windows, %d non-zero.",
                  nrow(wins), sum(wins[[signal_col]] > 0)))
  wins
}


# =============================================================================
# normalize_density()
# =============================================================================

normalize_density <- function(density_dt,
                              features  = NULL,
                              by_genome = FALSE) {
  
  dt   <- data.table::as.data.table(data.table::copy(density_dt))
  meta <- c("genome", "chr", "start", "end", "chr_len_bp", "other")
  
  if (is.null(features))
    features <- setdiff(names(dt), meta)
  
  if ("other" %in% features) {
    message("normalize_density: excluding 'other' (compositional remainder).")
    features <- setdiff(features, "other")
  }
  
  if (length(features) == 0L) {
    message("normalize_density: no features to normalise.")
    return(dt)
  }
  
  scope <- if (by_genome) "genome" else "chromosome"
  message(sprintf("normalize_density: minmax per %s on [%s]",
                  scope, paste(features, collapse = ", ")))
  
  for (feat in features) {
    if (!feat %in% names(dt)) {
      warning("Column not found, skipping: ", feat)
      next
    }
    
    rng <- range(dt[[feat]], na.rm = TRUE)
    message(sprintf("  '%s' raw: [%.4f, %.4f]", feat, rng[1], rng[2]))
    
    if (by_genome) {
      dt[, (feat) := {
        x  <- get(feat)
        mn <- min(x, na.rm = TRUE)
        mx <- max(x, na.rm = TRUE)
        if (is.na(mx) || mx == mn) rep(0, .N)
        else (x - mn) / (mx - mn)
      }, by = .(genome)]
    } else {
      dt[, (feat) := {
        x  <- get(feat)
        mn <- min(x, na.rm = TRUE)
        mx <- max(x, na.rm = TRUE)
        if (is.na(mx) || mx == mn) rep(0, .N)
        else (x - mn) / (mx - mn)
      }, by = .(genome, chr)]
    }
    
    rng_out <- range(dt[[feat]], na.rm = TRUE)
    message(sprintf("  '%s' normalised: [%.4f, %.4f]",
                    feat, rng_out[1], rng_out[2]))
  }
  
  dt
}


# =============================================================================
# prewindowed_to_density()
# =============================================================================

prewindowed_to_density <- function(
    dt,
    chr_sizes,
    value_col = "value",
    out_col   = NULL,
    normalize = TRUE,
    by_genome = FALSE
) {
  
  dt  <- data.table::as.data.table(data.table::copy(dt))
  cs  <- data.table::as.data.table(data.table::copy(chr_sizes))
  
  # ── validate ───────────────────────────────────────────────────────────────
  need <- c("genome", "chr", "start", "end")
  miss <- setdiff(need, names(dt))
  if (length(miss) > 0)
    stop("dt missing columns: ", paste(miss, collapse = ", "))
  if (!value_col %in% names(dt))
    stop(sprintf("value_col '%s' not found.\n  Available: %s",
                 value_col, paste(names(dt), collapse = ", ")))
  
  if (is.null(out_col)) out_col <- value_col
  
  # ── build result ───────────────────────────────────────────────────────────
  result <- dt[, .(genome,
                   chr,
                   start = as.numeric(start),
                   end   = as.numeric(end))]
  result[, (out_col) := as.numeric(dt[[value_col]])]
  
  # ── attach chr_len_bp ──────────────────────────────────────────────────────
  len_col <- .get_chr_len_col(cs)
  if (!"chr" %in% names(cs)) {
    alt <- intersect(c("chrom","chromosome","seqname","name"), names(cs))
    if (length(alt) == 0L) stop("chr_sizes has no 'chr' column.")
    data.table::setnames(cs, alt[1L], "chr")
  }
  cs_sub <- cs[, .(genome, chr, chr_len_bp = as.numeric(get(len_col)))]
  result  <- data.table::merge.data.table(result, cs_sub,
                                          by    = c("genome", "chr"),
                                          all.x = TRUE,
                                          sort  = FALSE)
  
  n_na <- sum(is.na(result$chr_len_bp))
  if (n_na > 0)
    warning(sprintf("%d windows have no chr_len_bp match.", n_na))
  
  # ── report ─────────────────────────────────────────────────────────────────
  message(sprintf(
    "prewindowed_to_density: %d windows, %d chr/genome combos",
    nrow(result),
    data.table::uniqueN(paste(result$genome, result$chr))))
  
  for (g in unique(result$genome)) {
    vals <- result[genome == g, get(out_col)]
    n_neg <- sum(vals < 0, na.rm = TRUE)
    message(sprintf(
      "  %s: raw [%.4f, %.4f] | %d non-zero | %d negative",
      g,
      min(vals, na.rm = TRUE),
      max(vals, na.rm = TRUE),
      sum(vals > 0, na.rm = TRUE),
      n_neg))
    if (n_neg > 0)
      message(sprintf(
        "    %d negative windows retained — will map to lower end of [0,1]",
        n_neg))
  }
  
  # ── normalize ──────────────────────────────────────────────────────────────
  if (!normalize) {
    message("prewindowed_to_density: normalize=FALSE, returning raw.")
    return(result)
  }
  
  normalize_density(result,
                    features  = out_col,
                    by_genome = by_genome)
}

# =============================================================================
# prepare_signal_density()
# =============================================================================

prepare_signal_density <- function(
    bedgraph,
    chr_sizes,
    genome_id  = NULL,
    win_size   = 1e5,
    step_size  = NULL,
    aggr       = c("mean", "sum", "max"),
    signal_col = "signal",
    normalize  = TRUE,
    by_genome  = FALSE
) {
  raw <- bedgraph_to_windows(
    bedgraph   = bedgraph,
    chr_sizes  = chr_sizes,
    genome_id  = genome_id,
    win_size   = win_size,
    step_size  = step_size,
    aggr       = match.arg(aggr),
    signal_col = signal_col
  )
  if (!normalize) return(raw)
  normalize_density(raw,
                    features  = signal_col,
                    by_genome = by_genome)
}

# =============================================================================
# prepare_feature_density()   BED → compute + normalize
# =============================================================================
prepare_feature_density <- function(
    features,
    chr_sizes,
    categories = c("gene", "LINE", "CRM", "other"),
    win_size   = 1e5,
    step_size  = NULL,
    method     = c("coverage", "count"),
    normalize  = TRUE,
    by_genome  = FALSE
) {
  raw <- calc_windowed_density(
    features   = features,
    chr_sizes  = chr_sizes,
    win_size   = win_size,
    step_size  = step_size,
    categories = categories,
    method     = match.arg(method)
  )
  if (!normalize) {
    message("prepare_feature_density: normalize=FALSE.")
    return(raw)
  }
  normalize_density(raw,
                    features  = setdiff(categories, "other"),
                    by_genome = by_genome)
}