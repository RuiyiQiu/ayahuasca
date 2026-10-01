library(ggplot2)
library(reshape2)
library(ggpubr)
library(tidyr)
library(dplyr)

# load data and process
df <- read.table("mapping_stats.txt", header=T, sep="\t", fill=T, na.strings = "NA")
df_long <- pivot_longer(df, cols = c(mean_depth, coverage, mapping_rate, SNP_number), 
                        names_to = "Metric", values_to = "Value")
# summarize basic statistics
summary(df)
## manually label on slide

# plot for overall metrics
ggplot(df_long, aes(x = Metric, y = Value, fill = Metric)) +
  geom_violin(alpha = 0.4) +
  geom_boxplot(width=0.2)+
  facet_wrap(~Metric, scales = "free", ncol=4) +
  theme_minimal() +
  theme(legend.position = "none",
        text = element_text(color = "black", size = 16),
        strip.text = element_text(size = 20, face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())
ggsave("mapping_stats.png", width = 12, height = 4)

# compare between 2 varieties
df_clean <- df %>%
  filter(!is.na(variety) & variety != "" & variety != " ")
target_metrics <- c("mapping_rate", "SNP_number")
df_compare <- pivot_longer(df_clean, cols = all_of(target_metrics), 
                           names_to = "Metric", values_to = "Value")

# summary
df_clean %>%
  group_by(variety) %>%
  summarise(mapping_rate=format(round(mean(mapping_rate, na.rm = TRUE),2)), # force print 2 digits
            SNP_number=mean(SNP_number, na.rm = TRUE))


ggplot(df_compare, aes(x = variety, y = Value, fill = variety)) +
  geom_violin(alpha = 0.4) +
  geom_boxplot(width=0.1)+
  scale_fill_manual(values = c("#EEE9BF","#A4D3EE"))+
  facet_wrap(~Metric, scales = "free", ncol=1) +
  stat_compare_means(method = "wilcox.test", 
                     label.y.npc = "top") +
  theme_minimal() +
  theme(legend.position = "none",
        text = element_text(color = "black", size = 16),
        strip.text = element_text(size = 18, face = "bold"),
        axis.title.x = element_blank())
ggsave("variety_stats.png", width = 4, height = 8)
