setwd("/Volumes/netscratch/dep_mercier/grp_marques/rqiu/Banisteropsis_Psychotria/Bcaapi_genome/TCN_v260520/plot")
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
chrlen_files <- c(TCN_h1 = "B_caapi_TCN_h1.chrlen",
                  TCN_h2 = "B_caapi_TCN_h2.chrlen",
                  CPR_h1 = "B_caapi_CPR_h1.chrlen",
                  CPR_h2 = "B_caapi_CPR_h2.chrlen")

chr_sizes  <- read_chrlen(chrlen_files)
gids <- c("TCN_h1", "TCN_h2", "CPR_h1", "CPR_h2")


##### 2.Reading the syri result #####
blk_syri <- syri_to_blkdt(
  syri_files = c("TCN_h1__TCN_h2" = "TCN_h1_TCN_h2syri.out",
                 "TCN_h2__CPR_h1" = "TCN_h2_CPR_h1syri.out",
                 "CPR_h1__CPR_h2" = "CPR_h1_CPR_h2syri.out"),
  genome_ids    = gids ,
  chr_sizes     = chr_sizes,
  minBlockSize = 10000)


##### 3.Reading genome annotation files #####
gene_files <- c(TCN_h1="TCN_h1.gene.bed", TCN_h2="TCN_h2.gene.bed", CPR_h1="CPR_h1.gene.bed", CPR_h2="CPR_h2.gene.bed")
te_files <- c(TCN_h1="TCN_h1.te.bed", TCN_h2="TCN_h2.te.bed", CPR_h1="CPR_h1.te.bed", CPR_h2="CPR_h2.te.bed")
crm_files <- c(TCN_h1="TCN_h1.crm.bed", TCN_h2="TCN_h2.crm.bed", CPR_h1="CPR_h1.crm.bed", CPR_h2="CPR_h2.crm.bed")
line_files   <- c(TCN_h1="TCN_h1.line.bed", TCN_h2="TCN_h2.line.bed", CPR_h1="CPR_h1.line.bed", CPR_h2="CPR_h2.line.bed")

genes_dt <- rbindlist(Map(read_bed, gene_files, names(gene_files), category = "bgene"))
te_dt <- rbindlist(Map(read_bed, te_files, names(te_files), category = "aTE"))
crm_dt <- rbindlist(Map(read_bed, crm_files, names(crm_files), category = "cCRM"))
line_dt   <- rbindlist(Map(read_bed, line_files, names(line_files), category = "LINE"))

all_features <- rbind(genes_dt, te_dt, crm_dt, line_dt)

cenh3_files <- c(TCN_h1="TCN_h1.cenh3.bedGraph", TCN_h2="TCN_h2.cenh3.bedGraph",
                 CPR_h1="CPR_h1.cenh3.bedGraph", CPR_h2="CPR_h2.cenh3.bedGraph")
cenh3_bg <- rbindlist(Map(read_bedgraph, cenh3_files, names(cenh3_files), category = "CENH3"))

cen_files <- c(TCN_h1 = "TCN_h1.cen.bed",
               TCN_h2 = "TCN_h2.cen.bed",
               CPR_h1 = "CPR_h1.cen.bed",
               CPR_h2 = "CPR_h2.cen.bed")
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
  normalize = FALSE)

track1 <- prepare_feature_density(
  features   = all_features,
  chr_sizes  = chr_sizes,
  categories = "LINE",
  win_size   = 1e6,
  method     = "coverage",
  normalize  = TRUE)

track2 <- prepare_feature_density(
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
  genomeIDs  = gids,
  showChrLabels = FALSE,
  space_above = 0,    # reserve space for density track (height 0.3 + buffer)
  space_below = 0.3,
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

##### Track 1 below: LINE #####
p2 <- add_feature_density_track(
  p,
  density_dt     = track1,
  features       = "LINE",
  feat_colors    = c(LINE = "#008585"),
  track_height   = 0.14,
  track_position = "below",
  plot_mode      = "overlay",
  track_gap      = 0,
  invert_y       = TRUE
)
p2

# Track 2 below: CRM (back) → other_TE → gene (front)
p3<- add_feature_density_track(
  p2,
  density_dt     = track2,
  features       = c("aTE", "bgene", "cCRM"),
  feat_colors    = c(aTE       = "#B4EEB4",
                     bgene     = "#1E90FF",
                     cCRM      = "#8B0000"),
  track_height   = 0.15,
  track_position = "below",
  plot_mode      = "overlay",
  track_gap      = 0.01,
  invert_y       = TRUE
)
p3


# CENH3 track above chromosomes
p4 <- add_feature_density_track(
  p3,
  density_dt     = cenh3,
  features       = "CENH3",
  feat_colors    = c(CENH3 = "#7B2D8B"),
  track_height   = 0.2,
  track_position = "above",
  plot_mode      = "overlay",
  invert_y       = FALSE)
p4


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

p5 <- p4 +
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
p5

ggsave(p5, file="plot.pdf", height=6, width=12)
