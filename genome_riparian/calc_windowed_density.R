# =============================================================================
# calc_windowed_density.R
#
# Robust, compatible windowed genomic feature density calculation.
# No non-equi joins, no type ambiguity, no name-collision bugs.
#
# Dependencies: data.table
# =============================================================================

calc_windowed_density <- function(
    features,          # data.table: genome, chr, start, end, category
    chr_sizes,         # data.table: genome, chr, size  (bp, from FAI)
    win_size   = 1e6,  # window size in bp
    step_size  = NULL, # step size in bp; NULL = non-overlapping (= win_size)
    categories = c("gene", "CRM", "TE", "other"),
    method     = "coverage"  # "coverage" (bp) or "count" (# features)
) {
  
  # ── 0. setup ───────────────────────────────────────────────────────────────
  if (!requireNamespace("data.table", quietly = TRUE))
    stop("data.table is required.")
  library(data.table)
  
  features  <- as.data.table(copy(features))
  chr_sizes <- as.data.table(copy(chr_sizes))
  
  # force all coordinates to plain numeric (avoids integer overflow & type mix)
  features[, start := as.numeric(start)]
  features[, end   := as.numeric(end)]
  chr_sizes[, size := as.numeric(size)]
  
  win_size  <- as.numeric(win_size)
  if (is.null(step_size)) step_size <- win_size
  step_size <- as.numeric(step_size)
  
  # "other" is derived, not looked up in features
  real_cats <- setdiff(categories, "other")
  do_other  <- "other" %in% categories
  
  # ── 1. validate inputs ─────────────────────────────────────────────────────
  req_feat <- c("genome", "chr", "start", "end", "category")
  miss <- setdiff(req_feat, names(features))
  if (length(miss) > 0)
    stop("features is missing columns: ", paste(miss, collapse = ", "))
  
  req_size <- c("genome", "chr", "size")
  miss2 <- setdiff(req_size, names(chr_sizes))
  if (length(miss2) > 0)
    stop("chr_sizes is missing columns: ", paste(miss2, collapse = ", "))
  
  # warn about categories not found
  found_cats <- unique(features$category)
  not_found  <- setdiff(real_cats, found_cats)
  if (length(not_found) > 0)
    message("Note: categories not found in features (will be 0): ",
            paste(not_found, collapse = ", "))
  
  # ── 2. build window table ──────────────────────────────────────────────────
  message("Building windows ...")
  
  win_list <- vector("list", nrow(chr_sizes))
  for (i in seq_len(nrow(chr_sizes))) {
    g    <- chr_sizes$genome[i]
    ch   <- chr_sizes$chr[i]
    sz   <- chr_sizes$size[i]
    
    ws   <- seq(1, sz, by = step_size)
    we   <- pmin(ws + win_size - 1, sz)
    
    win_list[[i]] <- data.table(
      genome          = g,
      chr             = ch,
      start           = ws,
      end             = we,
      win_size_actual = we - ws + 1   # true bp width (last window may be short)
    )
  }
  
  wins <- rbindlist(win_list)
  wins[, win_id := seq_len(.N)]
  
  # initialise all output columns to 0
  for (cat in categories)
    wins[, (cat) := 0]
  
  # ── 3. per-category, per-chromosome coverage sweep ─────────────────────────
  # Uses a simple two-pointer sweep after sorting — O(N log N), no joins.
  
  # helper: given sorted feature vectors (fs, fe) and a set of windows (ws, we),
  # return a numeric vector of covered bp per window
  .sweep_coverage <- function(fs, fe, ws, we) {
    nw  <- length(ws)
    nf  <- length(fs)
    out <- numeric(nw)       # plain numeric, no type issues
    if (nf == 0L) return(out)
    
    fi <- 1L                 # left pointer into features
    for (w in seq_len(nw)) {
      wstart <- ws[w]
      wend   <- we[w]
      
      # advance left pointer: skip features that end before this window
      while (fi <= nf && fe[fi] < wstart) fi <- fi + 1L
      if (fi > nf) break     # no more features for remaining windows
      
      # accumulate all features overlapping this window
      j     <- fi
      total <- 0
      while (j <= nf && fs[j] <= wend) {
        total <- total + pmin(fe[j], wend) - pmax(fs[j], wstart) + 1
        j <- j + 1L
      }
      out[w] <- total
    }
    out
  }
  
  # helper: count of features whose midpoint falls in each window
  .sweep_count <- function(fs, fe, ws, we) {
    mids <- (fs + fe) / 2
    nw   <- length(ws)
    out  <- numeric(nw)
    for (w in seq_len(nw)) {
      out[w] <- sum(mids >= ws[w] & mids <= we[w])
    }
    out
  }
  
  # unique (genome, chr) combinations present in wins
  chr_keys <- unique(wins[, .(genome, chr)])
  n_chr    <- nrow(chr_keys)
  
  for (cat in real_cats) {
    
    message(sprintf("  Processing category: %s", cat))
    
    f_cat <- features[category == cat,
                      .(genome, chr,
                        fs = as.numeric(start),
                        fe = as.numeric(end))]
    
    if (nrow(f_cat) == 0L) next
    
    for (ki in seq_len(n_chr)) {
      g  <- chr_keys$genome[ki]
      ch <- chr_keys$chr[ki]
      
      # windows for this chromosome
      widx <- wins[, genome == g & chr == ch]   # logical vector
      fw   <- wins[widx]
      if (nrow(fw) == 0L) next
      
      # features for this chromosome, sorted by start
      ff <- f_cat[genome == g & chr == ch]
      if (nrow(ff) == 0L) next
      ff <- ff[order(fs)]
      
      # compute coverage or count
      if (method == "coverage") {
        bp_vec <- .sweep_coverage(
          fs = ff$fs, fe = ff$fe,
          ws = fw$start, we = fw$end
        )
      } else {
        bp_vec <- .sweep_count(
          fs = ff$fs, fe = ff$fe,
          ws = fw$start, we = fw$end
        )
      }
      
      # write results back by win_id (safe, no row-order assumptions)
      target_ids <- fw$win_id
      wins[win_id %in% target_ids, (cat) := bp_vec]
    }
  }
  
  # ── 4. compute "other" as remainder ────────────────────────────────────────
  if (do_other) {
    message("  Computing 'other' as remainder ...")
    # sum of all real categories per window
    real_sums <- rowSums(wins[, real_cats, with = FALSE])
    wins[, other := pmax(0, win_size_actual - real_sums)]
  }
  
  # ── 5. clean up and return ─────────────────────────────────────────────────
  for (cat in categories) {
    wins[, (cat) := get(cat) / win_size_actual]
  }
  
  wins[, win_size_actual := NULL]
  wins[, win_id          := NULL]
  
  message("Done. ", nrow(wins), " windows across ",
          uniqueN(wins[, paste(genome, chr)]), " chromosomes.")
  wins
}
