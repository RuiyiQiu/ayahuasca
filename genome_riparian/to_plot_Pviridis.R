# R packages
library(data.table)
library(ggplot2)
#library(DEEPSPACE)

# functions
source("annotated_riparian.R")
source("read_files.R")
source("read_syri.R")
source("density_compute.R")
source("density_plot.R")


##### 1. Reading chr size #####
chrlen_files <- c(
  g1 = "g1.chrlen",
  g2 = "g2.chrlen",
  g3 = "g3.chrlen",
  g4 = "g4.chrlen",
  g5 = "g5.chrlen",
  g6 = "g6.chrlen")


chr_sizes  <- read_chrlen(chrlen_files)
gids <- c("g1", "g2", "g3", "g4", "g5", "g6")


##### 2.Reading the syri result #####
blk_syri <- syri_to_blkdt(
  syri_files = c("g1__g2" = "g1_g2syri.out",
                 "g2__g3" = "g2_g3syri.out",
                 "g3__g4" = "g3_g4syri.out",
                 "g4__g5" = "g4_g5syri.out",
                 "g5__g6" = "g5_g6syri.out"),
  genome_ids    = gids ,
  chr_sizes     = chr_sizes,
  minBlockSize = 1000)


##### 3.Reading genome annotation files #####
gene_files <- c(g1="g1.gene.bed", g2="g2.gene.bed", g3="g3.gene.bed", g4="g4.gene.bed", g5="g5.gene.bed", g6="g6.gene.bed")
crm_files <- c(g1="g1.crm.bed", g2="g2.crm.bed", g3="g3.crm.bed", g4="g4.crm.bed", g5="g5.crm.bed", g6="g6.crm.bed")
te_files   <- c(g1="g1.te.bed", g2="g2.te.bed", g3="g3.te.bed", g4="g4.te.bed", g5="g5.te.bed", g6="g6.te.bed") # TE except for CRM
genes_dt <- rbindlist(Map(read_bed, gene_files, names(gene_files), category = "bgene"))
crm_dt <- rbindlist(Map(read_bed, crm_files, names(crm_files), category = "cCRM"))
tes_dt   <- rbindlist(Map(read_bed, te_files,   names(te_files), category = "aTE"))

all_features <- rbind(genes_dt, crm_dt, tes_dt)

cenh3_files <- c(g1="g1.cenh3.bedGraph", g2="g2.cenh3.bedGraph", g3="g3.cenh3.bedGraph", g4="g4.cenh3.bedGraph", g5="g5.cenh3.bedGraph", g6="g6.cenh3.bedGraph")
cenh3_bg <- rbindlist(Map(read_bedgraph, cenh3_files, names(cenh3_files), category = "CENH3"))

cen_files <- c(g1 = "g1.cen.bed",
               g2 = "g2.cen.bed",
               g3 = "g3.cen.bed",
               g4 = "g4.cen.bed",
               g5 = "g5.cen.bed",
               g6 = "g6.cen.bed")
cen_dt <- rbindlist(Map(function(f, gid) {
  dt <- fread(f, header = FALSE, col.names = c("chr", "start", "end"))
  dt[, genome := gid]
  dt
}, cen_files, names(cen_files)))

##### 4. compute #####
cenh3 <- prewindowed_to_density(
  dt        = cenh3_bg,       # your table with genome|chr|start|end|value
  chr_sizes = chr_sizes,
  value_col = "value",
  out_col   = "CENH3", 
  normalize = TRUE)

track <- prepare_feature_density(
  features   = all_features,
  chr_sizes  = chr_sizes,
  categories = c("aTE", "bgene", "cCRM"),
  win_size   = 1e6,
  method     = "coverage",
  normalize  = TRUE)

##### 5. Plotting #####
p <- annotated_riparian(
  paf        = NULL,
  blk_dt     = blk_syri,
  minBlockSize        = 5e4,
  genomeIDs  = gids,
  showChrLabels = FALSE,
  space_above = 0,    # reserve space for density track (height 0.3 + buffer)
  space_below = 0.2,
  gapWidth   = 0.05,
  colorBy    = "type",
  typeColors = c(SYN   = "grey80",
                 INV   = "#FF7F00",
                 TRANS = "#00FF00",
                 DUP   = "#FF82AB"),
  braidAlpha = 0.5, 
  genomeLabelColor    = "black",
  bgColor             = "white")
p

# Get layout (already attached to p, or extract it)
layout <- attr(p, "chr_layout")
chrHeight <- attr(p, "chrHeight")

##### Centromeres #####
cen_layout <- merge(cen_dt,
                    layout[, c("genome", "chr", "xmin", "y")],
                    by = c("genome", "chr"), all.x = TRUE)

# xmin of the chromosome bar is the bp=0 origin for that chromosome
cen_layout[, `:=`(rect_xmin = xmin + start,
                  rect_xmax = xmin + end,
                  rect_ymin = y - chrHeight * 0.6,
                  rect_ymax = y + chrHeight * 0.6)]

# Add centromere layer on top of the existing plot
p <- p + geom_rect(data = as.data.frame(cen_layout),
                   mapping = aes(xmin = rect_xmin, xmax = rect_xmax,
                                 ymin = rect_ymin, ymax = rect_ymax),
                   fill    = "#E41A1C", colour  = NA, inherit.aes = FALSE)
p

# Track below
p2<- add_feature_density_track(
  p,
  density_dt     = track,
  features       = c("aTE", "bgene", "cCRM"),
  feat_colors    = c(aTE       = "#B4EEB4",
                     bgene     = "#1E90FF",
                     cCRM      = "#8B0000"),
  track_height   = 0.2,
  track_position = "below",
  plot_mode      = "overlay",
  track_gap      = 0.01,
  invert_y       = TRUE
)
p2


# CENH3 track above chromosomes
p3 <- add_feature_density_track(
  p2,
  density_dt     = cenh3,
  features       = "CENH3",
  feat_colors    = c(CENH3 = "#7B2D8B"),
  track_height   = 0.2,
  track_position = "above",
  plot_mode      = "overlay",
  invert_y       = FALSE)
p3

##### Scale bar #####
# Adjustable parameters
scalebar_mb    <- 50          # length in Mb
scalebar_color <- "black"     # line + label colour
scalebar_size  <- 0.5         # line thickness
scalebar_tsize <- 2.5         # label font size
scalebar_corner <- "topright"  # "bottomright", "bottomleft", "topright", "topleft"

layout     <- attr(p, "chr_layout")
n_genomes  <- length(unique(layout$genome))
bar_bp     <- scalebar_mb * 1e6

# x range of the plot (all chromosome bars)
x_range    <- range(c(layout$xmin, layout$xmax))
x_span     <- diff(x_range)

# y range: rows go from 0 to -(n_genomes-1), with 0.5 padding on each side
y_top      <-  0.5
y_bot      <- -(n_genomes - 1) - 0.5
y_span     <- y_top - y_bot

# Margin as fraction of axis range
x_margin   <- x_span * 0.02
y_margin   <- y_span * 0.05

# Anchor the bar
if (grepl("right",  scalebar_corner)) x_end   <- x_range[2] - x_margin  else
  x_end   <- x_range[1] + x_margin + bar_bp
x_start <- x_end - bar_bp

if (grepl("bottom", scalebar_corner)) y_bar   <- y_bot + y_margin  else
  y_bar   <- y_top - y_margin

label_text <- paste0(scalebar_mb, " Mb")

p4 <- p3 +
  # Horizontal bar
  annotate("segment",
           x = x_start, xend = x_end,
           y = y_bar,   yend = y_bar,
           colour    = scalebar_color,
           linewidth = scalebar_size) +
  # Left tick
  annotate("segment",
           x = x_start, xend = x_start,
           y = y_bar - y_span * 0.01, yend = y_bar + y_span * 0.01,
           colour    = scalebar_color,
           linewidth = scalebar_size) +
  # Right tick
  annotate("segment",
           x = x_end, xend = x_end,
           y = y_bar - y_span * 0.01, yend = y_bar + y_span * 0.01,
           colour    = scalebar_color,
           linewidth = scalebar_size) +
  # Label
  annotate("text",
           x      = (x_start + x_end) / 2,
           y      = y_bar + y_span * 0.025,
           label  = label_text,
           size   = scalebar_tsize,
           colour = scalebar_color,
           vjust  = 0)
p4

ggsave(p4, file="plot.pdf", height=6, width=12)
