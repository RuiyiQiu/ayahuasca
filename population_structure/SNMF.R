# install LEA R package
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("LEA")

library(LEA)
library(mapplots)
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(ggforce)
library(ggrepel)
my_colors <- c("#007AFF", "#FF3B30", "#34C759", "#AF52DE", "#FF9500", "#5AC8FA", "#FF72BE", "#FFCC00", "#4CD964")


##### read input #####
ped2geno("thinned_data.ped")

##### calculate #####
obj.at = snmf("thinned_data.geno", K = 1:10, ploidy = 2, entropy = T,
              CPU = 8, project = 'new')
#obj.at = load.snmfProject("thinned_data.snmfProject")
png(filename="cross-entropy.png")
plot(obj.at, col = "blue4", cex = 1.4, pch = 19)
dev.off()

##### Read Metadata #####
ids_raw <- read.table("id_order.txt", header = FALSE, stringsAsFactors = FALSE)
ids_raw$Library <- gsub("^.*/|\\.bam$", "", ids_raw$V1)
meta <- read.table('meta_2.txt', header = TRUE, sep = '\t')
ids_meta <- merge(ids_raw, meta, by = 'Library', all.x = TRUE, sort = FALSE)

##### Cumulative Matching Algorithm #####
ks_to_plot <- 2:8
max_k <- max(ks_to_plot)

# align the clustering between different Ks
aligned_q_list <- list()
for (k in ks_to_plot) {
  best_run = which.min(cross.entropy(obj.at, K = k))
  current_q <- as.data.frame(Q(obj.at, K = k, run = best_run))
  colnames(current_q) <- paste0("V", 1:k)
  
  if (k == min(ks_to_plot)) {
    # if only K=2, we just accept the default order
    colnames(current_q) <- paste0("C", 1:k)
    aligned_q_list[[as.character(k)]] <- current_q
  } else {
    # For subsequent Ks, we calculate which new column matches the previous K's columns
    prev_q <- aligned_q_list[[as.character(k-1)]]
    
    # Calculate Pearson correlation between all pairs of old and new clusters
    cor_matrix <- cor(current_q, prev_q)
    
    new_order <- numeric(k)
    used_prev <- numeric()
    # Match existing clusters based on highest correlation
    for(prev_idx in 1:(k-1)) {
      best_match <- which.max(cor_matrix[, prev_idx])
      new_order[prev_idx] <- best_match
      cor_matrix[best_match, ] <- -1 # Remove this cluster from further matching
    }
    
    # The newly appearing cluster in this K goes to the end
    new_order[k] <- setdiff(1:k, new_order[1:(k-1)])
    
    # Reorder the columns
    aligned_q <- current_q[, new_order]
    colnames(aligned_q) <- paste0("C", 1:k)
    aligned_q_list[[as.character(k)]] <- aligned_q
  }
}

##### barplot for multiple K ######
# sort by: 1. Dominant Cluster ID 2. Proportion
ref_q <- aligned_q_list[[as.character(max_k)]]
ref_q$Sample <- ids_meta$Sample
ref_q$dom_clust <- apply(ref_q[, 1:max_k], 1, which.max)
ref_q$max_val <- apply(ref_q[, 1:max_k], 1, max)
sample_ordering <- ref_q %>%
  arrange(dom_clust, max_val) %>%
  pull(Sample)

# merge oredered data for multiple K 
all_data <- data.frame()
for (k in ks_to_plot) {
  df_k <- aligned_q_list[[as.character(k)]]
  df_k$Sample <- ids_meta$Sample
  df_k$K_val <- paste0("K = ", k)
  
  df_long <- melt(df_k, id.vars = c("Sample", "K_val"), 
                  variable.name = "Cluster", value.name = "Admixture")
  all_data <- rbind(all_data, df_long)
}
all_data$Sample <- factor(all_data$Sample, levels = sample_ordering)
all_data$K_val <- factor(all_data$K_val, levels = paste0("K = ", ks_to_plot))

# plot
ggplot(all_data, aes(x = Sample, y = Admixture, fill = Cluster)) +
  geom_bar(stat = "identity", width = 1) +
  facet_grid(K_val ~ ., switch = "y") + 
  scale_fill_manual(values = my_colors) +
  theme_minimal() +
  labs(fill = "Ancestry") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        strip.text.y.left = element_text(angle = 0, face = "bold", size = 11),
        legend.position = "right",
        panel.grid = element_blank(),
        panel.spacing = unit(0.3, "lines"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3)) +
  scale_y_continuous(expand = c(0,0))
ggsave("Admixture_Kstacked_barplot.pdf", width = 20, height = 12)

##### Manually order by local and interior by ancestry #####
# sort by: 1. Local 2. Dominant Cluster ID 3. Proportion
my_local_order_vec <- read.table("local_order_2.txt")$V1
ref_q <- aligned_q_list[[as.character(max_k)]]
ref_q$Sample <- ids_meta$Sample
ref_q$Local <- ids_meta$Local
ref_q$dom_clust <- apply(ref_q[, 1:max_k], 1, which.max)
ref_q$max_val <- apply(ref_q[, 1:max_k], 1, max)
sample_ordering <- ref_q %>%
  mutate(Local = factor(Local, levels = my_local_order_vec)) %>%
  arrange(Local, dom_clust, max_val) %>%
  pull(Sample)

# merge oredered data for multiple K 
all_data <- data.frame()
for (k in ks_to_plot) {
  df_k <- aligned_q_list[[as.character(k)]]
  df_k$Sample <- ids_meta$Sample
  df_k$Local <- ids_meta$Local
  df_k$K_val <- paste0("K = ", k)
  
  df_long <- melt(df_k, id.vars = c("Sample", "Local", "K_val"), 
                  variable.name = "Cluster", value.name = "Admixture")
  all_data <- rbind(all_data, df_long)
}
all_data$Sample <- factor(all_data$Sample, levels = sample_ordering)
all_data$Local <- factor(all_data$Local, levels = my_local_order_vec)
all_data$K_val <- factor(all_data$K_val, levels = paste0("K = ", ks_to_plot))

# plotting
ggplot(all_data, aes(x = Sample, y = Admixture, fill = Cluster)) +
  geom_bar(stat = "identity", width = 1) +
  facet_grid(K_val ~ Local, switch = "y", space = "free_x", scales = "free_x") + 
  scale_fill_manual(values = my_colors) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        strip.text.y.left = element_text(angle = 0, face = "bold", size = 11),
        strip.text.x = element_text(size = 8, face = "bold"),
        panel.grid = element_blank(),
        panel.spacing = unit(0.1, "lines"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.1),
        legend.position = "right") +
  scale_y_continuous(expand = c(0,0))

ggsave("Admixture_Manual_Local_Ancestry.pdf", width = 20, height = 5)


## plotting 2
ggplot(all_data, aes(x = Sample, y = Admixture, fill = Cluster)) +
  geom_bar(stat = "identity", width = 1) +
  facet_grid(K_val ~ Local, switch = "y", space = "free_x", scales = "free_x") + 
  scale_fill_manual(values = my_colors) +
  theme_minimal(base_size = 14) +
  labs(fill = "Ancestry") +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        strip.text.y.left = element_text(angle = 0, face = "bold"),
        strip.text.x = element_text(angle = 90),
        panel.grid = element_blank(),
        panel.spacing = unit(0.1, "lines"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.1),
        legend.position = "right") +
  scale_y_continuous(expand = c(0,0))

ggsave("Admixture_Manual_Local_Ancestry.short.pdf", width = 12, height = 4)

###### choose K and plot on map ######
# K=5 or 8 seems interesting
best_K <- 5
best_K <- 8

# 1. Use the Aligned Q matrix from our previous steps
# This ensures "Cluster 1" in the pie chart is the same as "Cluster 1" in the barplot
q_best <- aligned_q_list[[as.character(best_K)]]

# 2. Integrate geographical information
# Use ids_meta because it already contains the cleaned Sample and Local info
data <- cbind(ids_meta, q_best)

# 3. Group by Local and calculate the mean ancestry for each location
pop_local <- data %>%
  group_by(Local) %>%
  summarise(
    Latitude = mean(Lat),   # Ensure column names match your meta file (Lat/Long)
    Longitude = mean(Long),
    # Select the Cluster columns we just added (C1, C2, etc.)
    across(starts_with("C"), mean) 
  )

# 4. Transform to long format for ggplot
pop_local_long <- melt(pop_local, id.vars = c("Local", "Latitude", "Longitude"),
                       variable.name = "Ancestry", value.name = "Coefficient")

# Unique coordinates for labels
pie_labels <- pop_local %>% distinct(Local, Longitude, Latitude)

# 5. Plot
ggplot() + 
  borders("world", colour="grey80", fill="grey95") + 
  # Adjust coordinates to fit your specific study area
  coord_fixed(xlim=c(min(data$Long)-2, max(data$Long)+2), 
              ylim=c(min(data$Lat)-2, max(data$Lat)+2)) +
  theme_bw() +
  # Pie charts
  geom_arc_bar(data = pop_local_long,
               aes(x0 = Longitude, y0 = Latitude, r0 = 0, r = 0.5, 
                   fill = Ancestry, amount = Coefficient),
               stat = "pie",
               inherit.aes = FALSE,
               color = "black", linewidth = 0.1) + 
  # Use the same color palette as the barplot
  scale_fill_manual(values = my_colors) +
  geom_text_repel(data = pie_labels, 
                  aes(x = Longitude, y = Latitude, label = Local),
                  size = 3,  
                  nudge_y = 0.5, force = 2) +
  labs(title = paste0("Population Admixture (K = ", best_K, ")"),
       fill = "Ancestry") +
  theme(panel.grid = element_blank(),
        text = element_text(size=12))

ggsave(paste0('Admixture_local_pie_map.K_', best_K, '.png'), height = 6, width = 7)


##### nudge the overlapping locals #####
pop_nudge <- pop_local_long %>%
  mutate(Plot_Long = Longitude, Plot_Lat = Latitude) %>%
  mutate(Plot_Long = ifelse(Local == "VP", Longitude + 1.5, Plot_Long),
         Plot_Lat  = ifelse(Local == "VP", Latitude, Plot_Lat)) %>%
  mutate(Plot_Long = ifelse(Local == "OP", Longitude, Plot_Long),
         Plot_Lat  = ifelse(Local == "OP", Latitude - 1.5, Plot_Lat)) %>%
  mutate(Plot_Long = ifelse(Local == "CS", Longitude + 1.5, Plot_Long),
         Plot_Lat  = ifelse(Local == "CS", Latitude - 0.5, Plot_Lat)) %>%
  mutate(Plot_Long = ifelse(Local == "TF", Longitude - 1.5, Plot_Long),
         Plot_Lat  = ifelse(Local == "TF", Latitude + 1.5, Plot_Lat)) %>%
  mutate(Plot_Long = ifelse(Local == "AF", Longitude, Plot_Long),
         Plot_Lat  = ifelse(Local == "AF", Latitude + 2, Plot_Lat)) %>%
  mutate(Plot_Long = ifelse(Local == "AF_2", Longitude + 1.5, Plot_Long),
         Plot_Lat  = ifelse(Local == "AF_2", Latitude - 1, Plot_Lat)) %>%
  mutate(Plot_Long = ifelse(Local == "TCN", Longitude - 1.5, Plot_Long),
         Plot_Lat  = ifelse(Local == "TCN", Latitude + 1.5, Plot_Lat)) %>%
  mutate(Plot_Long = ifelse(Local == "TCN_ref", Longitude + 1.5, Plot_Long),
         Plot_Lat  = ifelse(Local == "TCN_ref", Latitude + 1.5, Plot_Lat)) %>%
  mutate(Plot_Long = ifelse(Local == "CPR", Longitude - 1.5, Plot_Long),
         Plot_Lat  = ifelse(Local == "CPR", Latitude - 1.5, Plot_Lat)) %>%
  mutate(Plot_Long = ifelse(Local == "CPR_ref", Longitude + 1.5, Plot_Long),
         Plot_Lat  = ifelse(Local == "CPR_ref", Latitude - 1.5, Plot_Lat)) %>%
  mutate(Plot_Long = ifelse(Local == "RA", Longitude - 1.5, Plot_Long),
         Plot_Lat  = ifelse(Local == "RA", Latitude + 1.5, Plot_Lat)) %>%
  mutate(Plot_Long = ifelse(Local == "SL", Longitude + 0.5, Plot_Long),
         Plot_Lat  = ifelse(Local == "SL", Latitude + 1.5, Plot_Lat)) %>%
  mutate(Plot_Long = ifelse(Local == "PC", Longitude - 1.5, Plot_Long),
         Plot_Lat  = ifelse(Local == "PC", Latitude + 1.5, Plot_Lat)) %>%
  mutate(Plot_Long = ifelse(Local == "NE", Longitude + 2, Plot_Long),
         Plot_Lat  = ifelse(Local == "NE", Latitude - 0.5, Plot_Lat)) %>%
  mutate(Plot_Long = ifelse(Local == "EX", Longitude, Plot_Long),
         Plot_Lat  = ifelse(Local == "EX", Latitude + 1.5, Plot_Lat)) %>%
  mutate(Plot_Long = ifelse(Local == "AC", Longitude, Plot_Long),
         Plot_Lat  = ifelse(Local == "AC", Latitude - 1.5, Plot_Lat)) %>%
  mutate(Plot_Long = ifelse(Local == "CX", Longitude - 1.5, Plot_Long),
         Plot_Lat  = ifelse(Local == "CX", Latitude - 1.5, Plot_Lat))%>%
  mutate(Plot_Long = ifelse(Local == "TR", Longitude - 1.5, Plot_Long),
         Plot_Lat  = ifelse(Local == "TR", Latitude - 1.5, Plot_Lat))

ggplot() + 
  borders("world", colour="grey80", fill="grey95") + 
  coord_fixed(xlim=c(min(data$Long)-2, max(data$Long)+5), 
              ylim=c(min(data$Lat)-2, max(data$Lat)+2)) +
  theme_bw(base_size = 14) +
  
  # 2. Add Leader Lines (Segments)
  # This draws a thin line from the true location to the center of the shifted pie
  geom_segment(data = pop_nudge %>% distinct(Local, .keep_all = TRUE),
               aes(x = Plot_Long, y = Plot_Lat, 
                   xend = Longitude, yend = Latitude),
               colour = "grey40", 
               linewidth = 0.4,
               arrow = arrow(length = unit(0.1, "cm"), type = "closed"),
               alpha = 0.7) +
  
  # 3. Pie charts (now using Plot_Long/Lat)
  geom_arc_bar(data = pop_nudge,
               aes(x0 = Plot_Long, y0 = Plot_Lat, r0 = 0, r = 1, 
                   fill = Ancestry, amount = Coefficient),
               stat = "pie",
               inherit.aes = FALSE,
               color = "black", linewidth = 0.1) + 
  
  scale_fill_manual(values = my_colors) +
  
  # 4. Repel Labels (linked to the shifted Pie position)
  geom_text_repel(data = pop_nudge %>% distinct(Local, .keep_all = TRUE), 
                  aes(x = Plot_Long, y = Plot_Lat, label = Local),
                  size = 3, 
                  force = 1) +
  
  labs(fill = "Ancestry") +
  xlab("Longitude") +
  ylab("Latitude")+
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "none")

ggsave(paste0('Admixture_local_pie_map.nudged.K_', best_K, '.pdf'), height = 6, width = 10)
