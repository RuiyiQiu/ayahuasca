# =============================================================================
# density_plot.R
#
# add_feature_density_track()
#   Universal track function. Handles both overlay and stack modes,
#   both BED-derived density and bedGraph-derived signal.
#
# Track placement is automatic via the registry system.
# =============================================================================

utils::globalVariables(c(
  "x_mid","y_base","group_id","fill_col","feat",
  "reg_start","reg_end","start_clipped","end_clipped","xmin","xmax","chr_len"
))

`%||%` <- function(a, b) if (!is.null(a)) a else b


# =============================================================================
# .compute_feature_areas()
# =============================================================================

.compute_feature_areas <- function(density_dt,
                                   chr_layout,
                                   features,
                                   feat_colors,
                                   track_height,
                                   y_base_dt,
                                   direction,
                                   invert_y,
                                   plot_mode,
                                   smooth_method,
                                   smooth_n,
                                   smooth_span) {
  
  lay    <- data.table::as.data.table(chr_layout)
  dt     <- data.table::as.data.table(data.table::copy(density_dt))
  yb     <- data.table::as.data.table(y_base_dt)
  
  # ── direction + invert ─────────────────────────────────────────────────────
  #
  #  direction  invert_y   meaning
  #  ─────────  ────────   ──────────────────────────────────────────────────
  #  above      FALSE      baseline at chr top,    peak grows UP   (normal)
  #  above      TRUE       baseline at outer top,  peak grows DOWN toward chr
  #  below      FALSE      baseline at chr bottom, peak grows DOWN (normal)
  #  below      TRUE       baseline at outer bottom, peak grows UP toward chr
  #
  #  invert_y=TRUE shifts y_base to the OUTER edge of the track,
  #  and flips sign so the signal grows back toward the chromosome.
  
  natural_sign <- if (direction == "above") +1L else -1L
  
  if (invert_y) {
    # baseline is outer edge: y_base + natural_sign * track_height
    # signal grows in opposite direction (toward chr)
    sign_d        <- -natural_sign
    base_offset   <- natural_sign * track_height   # shift baseline to outer edge
  } else {
    # baseline is inner edge (chr edge): no shift
    sign_d        <- natural_sign
    base_offset   <- 0
  }
  
  # ── layout join ────────────────────────────────────────────────────────────
  has_reg <- "reg_start" %in% names(lay)
  if (has_reg) {
    lay_sub <- lay[, .(genome, chr, xmin, xmax,
                       lay_chr_len = chr_len,
                       reg_start, reg_end)]
  } else {
    lay_sub <- lay[, .(genome, chr, xmin, xmax,
                       lay_chr_len = chr_len)]
    lay_sub[, `:=`(reg_start = 0L, reg_end = lay_chr_len)]
  }
  
  dt <- data.table::merge.data.table(dt, lay_sub,
                                     by = c("genome","chr"), all.x = TRUE)
  dt <- data.table::merge.data.table(dt, yb,
                                     by = c("genome","chr"), all.x = TRUE)
  
  dt <- dt[!is.na(xmin) & !is.na(y_base)]
  dt <- dt[end >= reg_start & start <= reg_end]
  if (nrow(dt) == 0L) {
    warning(".compute_feature_areas: no windows in visible region.")
    return(.empty_area_dt())
  }
  
  dt[, start_clipped := pmax(start, reg_start)]
  dt[, end_clipped   := pmin(end,   reg_end)]
  dt[, x_mid := xmin + ((start_clipped + end_clipped) / 2 - reg_start)]
  dt <- dt[x_mid >= xmin & x_mid <= xmax]
  if (nrow(dt) == 0L) {
    warning(".compute_feature_areas: no windows after x-range filter.")
    return(.empty_area_dt())
  }
  
  # apply base_offset: shift y_base to outer edge when invert_y=TRUE
  dt[, y_base := y_base + base_offset]
  
  # ── smoother ───────────────────────────────────────────────────────────────
  .smooth <- function(x_raw, y_raw, x_out, method, span) {
    if (length(unique(y_raw)) == 1L)
      return(rep(y_raw[1L], length(x_out)))
    if (method == "spline") {
      tryCatch(
        stats::spline(x_raw, y_raw, xout = x_out, method = "natural")$y,
        error = function(e)
          stats::approx(x_raw, y_raw, xout = x_out, rule = 2)$y)
    } else {
      tryCatch({
        fit <- stats::loess(y ~ x,
                            data    = data.frame(x = x_raw, y = y_raw),
                            span    = span, degree = 2,
                            control = stats::loess.control(surface = "direct"))
        stats::predict(fit, newdata = data.frame(x = x_out))
      }, error = function(e)
        stats::approx(x_raw, y_raw, xout = x_out, rule = 2)$y)
    }
  }
  
  poly_list     <- list()
  group_counter <- 0L
  chr_groups    <- unique(dt[, .(genome, chr, xmin, xmax, y_base)])
  
  # ── OVERLAY ────────────────────────────────────────────────────────────────
  if (plot_mode == "overlay") {
    
    for (feat in features) {
      if (!feat %in% names(dt)) next
      
      # safe color lookup — use first available color as fallback
      col_fill  <- if (!is.null(feat_colors) && feat %in% names(feat_colors))
        feat_colors[[feat]]
      else if (!is.null(feat_colors) && length(feat_colors) > 0)
        feat_colors[[1]]
      else "grey60"
      
      for (ci in seq_len(nrow(chr_groups))) {
        cg     <- chr_groups[ci]
        dt_chr <- dt[genome == cg$genome & chr == cg$chr]
        if (nrow(dt_chr) < 2L) next
        data.table::setorder(dt_chr, x_mid)
        
        x_raw <- dt_chr$x_mid
        y_raw <- dt_chr[[feat]]
        x_out <- seq(min(x_raw), max(x_raw), length.out = smooth_n)
        
        y_sm        <- pmax(0, pmin(1, .smooth(x_raw, y_raw, x_out,
                                               smooth_method, smooth_span)))
        y_signal    <- cg$y_base + sign_d * y_sm * track_height
        y_base_line <- rep(cg$y_base, length(x_out))
        
        group_counter <- group_counter + 1L
        gid <- paste(cg$genome, cg$chr, feat, group_counter, sep = "__")
        poly_list[[gid]] <- data.table::data.table(
          x          = c(x_out,       rev(x_out)),
          y          = c(y_signal,     rev(y_base_line)),
          fill_col   = col_fill,
          feat       = feat,
          group_id   = gid
        )
      }
    }
    
    # ── STACK ──────────────────────────────────────────────────────────────────
  } else {
    
    for (ci in seq_len(nrow(chr_groups))) {
      cg     <- chr_groups[ci]
      dt_chr <- dt[genome == cg$genome & chr == cg$chr]
      if (nrow(dt_chr) < 2L) next
      data.table::setorder(dt_chr, x_mid)
      
      x_raw  <- dt_chr$x_mid
      x_out  <- seq(min(x_raw), max(x_raw), length.out = smooth_n)
      cum_sm <- rep(0, length(x_out))
      
      for (feat in features) {
        if (!feat %in% names(dt_chr)) next
        
        col_fill  <- if (!is.null(feat_colors) && feat %in% names(feat_colors))
          feat_colors[[feat]]
        else if (!is.null(feat_colors) && length(feat_colors) > 0)
          feat_colors[[1]]
        else "grey60"
        
        y_raw    <- dt_chr[[feat]]
        y_sm     <- pmax(0, .smooth(x_raw, y_raw, x_out,
                                    smooth_method, smooth_span))
        y_bottom <- cg$y_base + sign_d * cum_sm          * track_height
        y_top    <- cg$y_base + sign_d * (cum_sm + y_sm) * track_height
        
        group_counter <- group_counter + 1L
        gid <- paste(cg$genome, cg$chr, feat, group_counter, sep = "__")
        poly_list[[gid]] <- data.table::data.table(
          x          = c(x_out,    rev(x_out)),
          y          = c(y_top,     rev(y_bottom)),
          fill_col   = col_fill,
          feat       = feat,
          group_id   = gid
        )
        cum_sm <- cum_sm + y_sm
      }
    }
  }
  
  if (length(poly_list) == 0L) return(.empty_area_dt())
  data.table::rbindlist(poly_list)
}

.empty_area_dt <- function()
  data.table::data.table(x=numeric(), y=numeric(), fill_col=character(),
                         feat=character(), fill_alpha=numeric(),
                         group_id=character())


# =============================================================================
# add_feature_density_track()
# =============================================================================

add_feature_density_track <- function(
    p,
    density_dt,
    features,
    feat_colors     = NULL,
    track_height    = 0.3,
    track_position  = c("below", "above"),
    invert_y        = FALSE,
    plot_mode       = c("overlay", "stack"),
    smooth_method   = c("spline", "loess"),
    smooth_n        = 500,
    smooth_span     = 0.15,
    track_gap       = 0,
    add_baseline    = TRUE,
    baseline_color  = "grey40",
    baseline_size   = 0.3,
    outline_color   = NA,
    outline_size    = 0.2
) {
  
  track_position <- match.arg(track_position)
  plot_mode      <- match.arg(plot_mode)
  smooth_method  <- match.arg(smooth_method)
  direction      <- track_position
  
  chr_layout  <- attr(p, "chr_layout")
  chr_height  <- attr(p, "chrHeight")      %||% 0.05
  y_buffer    <- attr(p, "yBuffer")        %||% 0.02
  space_above <- attr(p, "space_above")    %||% 0
  space_below <- attr(p, "space_below")    %||% 0
  registry    <- attr(p, "track_registry") %||% list()
  
  if (is.null(chr_layout))
    stop("No chr_layout found. p must come from annotated_riparian().")
  
  dt   <- data.table::as.data.table(density_dt)
  miss <- setdiff(features, names(dt))
  if (length(miss) > 0)
    stop("Features not in density_dt: ", paste(miss, collapse = ", "),
         "\n  Available: ", paste(names(dt), collapse = ", "))
  if (!"chr_len_bp" %in% names(dt))
    stop("density_dt needs 'chr_len_bp'. Use prewindowed_to_density() ",
         "or prepare_feature_density().")
  
  # ── fix color names if mismatch ─────────────────────────────────────
  if (!is.null(feat_colors) &&
      (is.null(names(feat_colors)) ||
       !any(features %in% names(feat_colors)))) {
    message("  Assigning feat_colors by position.")
    names(feat_colors) <- features[seq_along(feat_colors)]
  }
  
  # ── defaults ───────────────────────────────────────────────────────────────
  default_pal <- c("#4393C3","#D6604D","#F4A582","#66C2A5",
                   "#FC8D62","#8DA0CB","#BEBEBE","#A6D854")
  if (is.null(feat_colors))
    feat_colors <- stats::setNames(
      default_pal[seq_along(features)], features)
  else
    stats::setNames(rep(0.75, length(features)), features)
  
  # ── report data range (informational only, no warning) ────────────────────
  for (feat in features) {
    rng <- range(dt[[feat]], na.rm = TRUE)
    message(sprintf(
      "  '%s': range [%.4f, %.4f]%s",
      feat, rng[1], rng[2],
      if (rng[2] > 1.001) " (raw, will be scaled to track height)" else ""))
  }
  
  if (plot_mode == "stack") {
    rs <- rowSums(as.matrix(dt[, features, with = FALSE]), na.rm = TRUE)
    if (mean(rs <= 1.01, na.rm = TRUE) < 0.9)
      warning("plot_mode='stack': many windows sum > 1. ",
              "Stack expects raw proportions.")
  }
  
  message(sprintf(
    "add_feature_density_track: [%s] | %s | mode=%s | invert_y=%s | h=%.2f",
    paste(features, collapse="+"), direction,
    plot_mode, invert_y, track_height))
  
  y_base_dt <- .get_y_base(
    chr_layout = chr_layout, chr_height = chr_height,
    y_buffer   = y_buffer,   direction  = direction,
    registry   = registry,   track_gap  = track_gap)
  
  areas <- .compute_feature_areas(
    density_dt    = dt,      chr_layout    = chr_layout,
    features      = features, feat_colors  = feat_colors,
    track_height = track_height,
    y_base_dt     = y_base_dt, direction    = direction,
    invert_y      = invert_y,  plot_mode    = plot_mode,
    smooth_method = smooth_method, smooth_n = smooth_n,
    smooth_span   = smooth_span)
  
  if (nrow(areas) == 0L) {
    warning("add_feature_density_track: no polygons produced.")
    return(p)
  }
  message(sprintf("  %d points, %d groups",
                  nrow(areas), data.table::uniqueN(areas$group_id)))
  
  if (add_baseline) {
    lay <- data.table::as.data.table(chr_layout)
    bl  <- data.table::merge.data.table(
      lay[, .(genome, chr, xmin, xmax)], y_base_dt, by = c("genome","chr"))
    natural_sign <- if (direction == "above") +1L else -1L
    base_offset  <- if (invert_y) natural_sign * track_height else 0
    bl[, y_base  := y_base + base_offset]
    p <- p + ggplot2::geom_segment(
      data        = as.data.frame(bl),
      mapping     = ggplot2::aes(x=xmin, xend=xmax, y=y_base, yend=y_base),
      colour      = baseline_color, linewidth = baseline_size,
      inherit.aes = FALSE)
  }
  
  for (feat in features) {
    fa <- areas[areas$feat == feat, ]
    if (nrow(fa) == 0L) next
    p  <- p + ggplot2::geom_polygon(
      data        = as.data.frame(fa),
      mapping     = ggplot2::aes(x=x, y=y, group=group_id, fill=fill_col),
      colour = outline_color, linewidth = outline_size,
      inherit.aes = FALSE)
  }
  p <- p + ggplot2::scale_fill_identity()
  
  registry[[length(registry) + 1L]] <- list(
    direction = direction, track_height = track_height, track_gap = track_gap)
  attr(p, "chr_layout")     <- chr_layout
  attr(p, "chrHeight")      <- chr_height
  attr(p, "yBuffer")        <- y_buffer
  attr(p, "space_above")    <- space_above
  attr(p, "space_below")    <- space_below
  attr(p, "track_registry") <- registry
  p
}
