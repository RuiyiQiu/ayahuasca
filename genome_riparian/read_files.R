read_gff3 <- function(gff_file, genome_id,
                      feature_types = c("gene", "transposable_element",
                                        "repeat_region", "mobile_genetic_element")) {
  # skip comment lines
  gff <- fread(gff_file, sep = "\t", header = FALSE, skip = "#",
               col.names = c("chr","source","type","start","end",
                             "score","strand","phase","attributes"))
  
  # keep only desired feature types (case-insensitive)
  gff <- gff[tolower(type) %in% tolower(feature_types)]
  
  # broad category
  gff[, category := fcase(
    tolower(type) %in% c("gene","mrna","pseudogene"),             "gene",
    tolower(type) %in% c("transposable_element",
                         "mobile_genetic_element",
                         "repeat_region",
                         "retrotransposon",
                         "dna_transposon"),                       "TE",
    default = "unknown"
  )]
  
  gff[, .(genome = genome_id, chr, start, end, category)]
}

read_bed <- function(bed_file, genome_id, category = "gene") {
  # Read without col.names first to count columns dynamically
  bed <- fread(bed_file, header = FALSE)
  
  # Name only the columns we care about; extras get generic names V4, V5...
  ncols <- ncol(bed)
  base_names <- c("chr", "start", "end", "name", "score", "strand",
                  "source", "type", "phase", "attributes")
  new_names <- base_names[seq_len(min(ncols, length(base_names)))]
  if (ncols > length(base_names))
    new_names <- c(new_names, paste0("V", seq(length(base_names)+1, ncols)))
  setnames(bed, new_names)
  
  # BED is 0-based half-open → convert to 1-based closed
  bed[, start := start + 1L]
  bed[, `:=`(genome = genome_id, category = category)]
  bed[, .(genome, chr, start, end, category)]
}

read_bedgraph <- function(bg_file, genome_id, category = "CENH3") {
  bg <- fread(bg_file, header = FALSE,
              col.names = c("chr", "start", "end", "value"))  # name it directly here
  bg[, `:=`(start  = start + 1L,
            genome = genome_id)]
  bg[, .(genome, chr, start, end, value, category)]   # select by name safely
}

read_fai <- function(fai_file, genome_id) {
  fai <- fread(fai_file, header = FALSE,
               col.names = c("chr","size","offset","bases_per_line","bytes_per_line"))
  fai[, .(genome = genome_id, chr, size)]
}

read_chromsizes <- function(f, genome_id) {
  dt <- fread(f, header = FALSE, col.names = c("chr","size"))
  dt[, genome := genome_id]
  dt[, .(genome, chr, size)]
}

