calc_windowed_density_fast <- function(features,
                                       chr_sizes,
                                       win_size   = 1e6,
                                       step_size  = NULL,
                                       categories = c("gene", "CRM", "TE", "other"),
                                       method     = "coverage") {
  
  features  <- as.data.table(features)
  chr_sizes <- as.data.table(chr_sizes)
  if (is.null(step_size)) step_size <- win_size
  
  # categories to actually query from features (exclude "other" — computed at end)
  real_cats <- setdiff(categories, "other")
  do_other  <- "other" %in% categories
  
  # ── 1. build window table ──────────────────────────────────────────────────
  wins <- chr_sizes[, {
    ws <- seq(1L, size, by = step_size)
    we <- pmin(ws + as.integer(win_size) - 1L, size)
    .(start = ws, end = we, win_size_actual = we - ws + 1L)
  }, by = .(genome, chr)]
  wins[, win_id := .I]
  
  # initialise all category columns to 0
  for (cat in categories) wins[, (cat) := 0L]
  
  # ── 2. for each real category, compute per-window coverage with a simple
  #        vectorised approach (no non-equi join — avoids the "each" bug) ────
  
  for (cat in real_cats) {
    
    f_cat <- features[category == cat,
                      .(genome, chr,
                        fs = as.integer(start),
                        fe = as.integer(end))]
    if (nrow(f_cat) == 0) next
    
    # process per (genome, chr) chunk — avoids any cross-chr contamination
    chr_chunks <- unique(wins[, .(genome, chr)])
    
    for (row_i in seq_len(nrow(chr_chunks))) {
      g  <- chr_chunks$genome[row_i]
      ch <- chr_chunks$chr[row_i]
      
      fw <- wins[genome == g & chr == ch]          # windows for this chr
      ff <- f_cat[genome == g & chr == ch]         # features for this chr
      
      if (nrow(ff) == 0) next
      
      # vectorised overlap: outer comparison matrix would be too large for
      # big chromosomes, so we use an interval sweep instead
      
      # sort features by start
      setorder(ff, fs)
      
      # for each window, sum intersected bp using findInterval
      fs_vec <- ff$fs
      fe_vec <- ff$fe
      
      bp_per_win <- vapply(seq_len(nrow(fw)), function(w) {
        ws <- fw$start[w]
        we <- fw$end[w]
        # features that overlap [ws, we]: fs <= we AND fe >= ws
        idx <- which(fs_vec <= we & fe_vec >= ws)
        if (length(idx) == 0L) return(0L)
        sum(pmin(fe_vec[idx], we) - pmax(fs_vec[idx], ws) + 1L)
      }, integer(1))
      
      # write back into wins
      wins[genome == g & chr == ch, (cat) := bp_per_win]
    }
  }
  
  # ── 3. compute "other" as remainder ───────────────────────────────────────
  if (do_other) {
    wins[, other := pmax(0L, win_size_actual - rowSums(.SD)),
         .SDcols = real_cats]
  }
  
  wins[, win_size_actual := NULL]
  wins[, win_id          := NULL]
  wins
}
