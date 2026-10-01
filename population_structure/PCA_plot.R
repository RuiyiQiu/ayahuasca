library(ggplot2)
library(ggrepel)
library(scatterplot3d)
library(scales)
prefix <- "merge_bi_SNP.flt.2"

# Read the eigenval file
eigenval_file <- paste0(prefix, ".eigenval")
eigenval <- read.table(eigenval_file, header = F, sep = "\t", comment.char = "#")
colnames(eigenval) <- c("eval", "eval_percent")

# Read the eigenvec file
pca_file <- paste0(prefix, ".eigenvec")
PCA <- read.table(pca_file, header = TRUE)
PCA$Sample <- gsub("^.*/|\\.bam$", "", PCA$SampleName)  # Clean sample names

# Read local information
local <- read.table('sample_local.all.txt',header = T, sep = '\t')
PCA_local <- merge(local, PCA, by.x='Library', by.y='Sample')

PCA_local$Local <- as.factor(PCA_local$Local)
PCA_local$Cluster <- as.factor(PCA_local$Cluster)

# process colors and shapes
unique_locals <- levels(PCA_local$Local)
color_palette <- scales::hue_pal()(length(unique_locals))
names(color_palette) <- unique_locals

transparent_palette <- scales::alpha(color_palette, 0.6)

print(data.frame(Local = names(color_palette), Color = color_palette, row.names = NULL))
write.table(data.frame(Local = names(color_palette), Color = color_palette), 
            file = "Local_Color_Mapping.txt", sep = "\t", row.names = FALSE, quote = FALSE)

unique_clusters <- levels(PCA_local$Cluster)
shape_values <- c(21, 22, 23, 24, 25)[1:length(unique_clusters)]
names(shape_values) <- unique_clusters

pt_colors <- transparent_palette[as.character(PCA_local$Local)]
pt_shapes <- shape_values[as.character(PCA_local$Cluster)]

# axis label
x_label <- paste0("PC1 (", round(eigenval$eval_percent[1], 2), "%)")
y_label <- paste0("PC2 (", round(eigenval$eval_percent[2], 2), "%)")
z_label <- paste0("PC3 (", round(eigenval$eval_percent[3], 2), "%)")


# plot 3D PCA
pdf("PCA_3D_Vector_MainPlot.pdf", width = 8, height = 8)
s3d <- scatterplot3d(x = PCA_local$PC1, y = PCA_local$PC2, z = PCA_local$PC3,
                     xlab = x_label, ylab = y_label, zlab = z_label,
                     color = pt_colors, bg = pt_colors, pch = pt_shapes, cex.symbols = 1.5,
                     grid = TRUE, box = TRUE, angle = 55)
dev.off()

# plot legeneds
pdf("PCA_3D_Legends_Only.pdf", width = 8, height = 8)
plot.new()
legend("topleft", 
       legend = names(transparent_palette),
       fill = transparent_palette, 
       border = transparent_palette,
       title = "Local", 
       cex = 1.2, 
       bty = "n")
legend("topright", 
       legend = names(shape_values),
       pch = shape_values, 
       col = scales::alpha("gray50", 0.6),
       pt.bg = scales::alpha("gray50", 0.6),
       title = "Cluster", 
       cex = 1.2, 
       bty = "n")
dev.off()




# === Specify PCs to Plot ===
x_PC <- "PC2"  # Replace with the desired x-axis PC
y_PC <- "PC3"  # Replace with the desired y-axis PC

# Extract Variance Explained% for the selected PCs
x_PC_index <- as.numeric(gsub("PC", "", x_PC))  # Get PC index (e.g., PC1 -> 1)
y_PC_index <- as.numeric(gsub("PC", "", y_PC))

x_var_explained <- eigenval$eval_percent[x_PC_index]
y_var_explained <- eigenval$eval_percent[y_PC_index]

# === Create the Plot ===
# Dynamically adjust labels for axes
x_label <- paste0(x_PC, " (", round(x_var_explained, 2), "%)")
y_label <- paste0(y_PC, " (", round(y_var_explained, 2), "%)")

# refine PCA plotting
ggplot(data = PCA_local, mapping = aes_string(x = x_PC, y = y_PC, 
                                            color = "Local", fill = "Local",
                                            shape = "Cluster")) +
  geom_point(size = 3)+
  labs(x = x_label, y = y_label) +
  scale_shape_manual(values = c(21, 22, 23, 24, 25)) +
  scale_color_manual(values = transparent_palette)+
  scale_fill_manual(values = transparent_palette)+
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        legend.position = "none")
# Save the plot with dynamic filename
ggsave(paste0(x_PC, "_", y_PC, ".cluster_local", ".pdf"), height=5, width=5)

# ---- Plot Distribution by PCA cluster ----
library(raster)
library(terra)
library(viridis)
ggplot() + 
  borders("world", colour="grey80", fill="grey95") + 
  coord_fixed(xlim=c(min(PCA_local$Long)-2, max(PCA_local$Long)+2), 
              ylim=c(min(PCA_local$Lat)-2, max(PCA_local$Lat)+2)) +
  geom_point(data=PCA_local,
             aes(x=Long,y=Lat, 
                 fill=Local, color=Local,
                 shape=Cluster), size = 3)+
  scale_color_manual(values = transparent_palette)+
  scale_fill_manual(values = transparent_palette)+
  scale_shape_manual(values = c(21, 22, 23, 24, 25)) + #same as in PCA plot
  theme_bw(base_size = 14) +
  xlab("Longitude") +
  ylab("Latitude")+
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "none")
ggsave(file='PCAcluster_distribution.pdf', height = 6, width = 10)

