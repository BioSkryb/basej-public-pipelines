library(ggplot2)
library(scales)
library(dplyr)

# Get the inputs
args <- commandArgs(trailingOnly=TRUE)
file_name <- args[1]

#Read input
ado <- read.table(file_name,header = T,sep = "\t")

# Plot the interval summary
plot <- ado %>%
  mutate(interval = sub("^.*_([^_]+)$", "\\1", File_Interval)) %>%
  ggplot(.,aes(interval,Prop)) +
  geom_boxplot()+
  labs(
    y="% Allelic Balance",
    title="% Allelic Balance Comparison",
    fill="Coverage",
    x = "Coverage"
  ) +
  theme_minimal()+
  theme(
    plot.title = element_text(size = 8),
    strip.text.x = element_text(size=4),
    axis.text.x = element_text(size = 4, angle = 90),
    axis.text.y = element_text(size = 4),
    axis.title.x = element_text(size = 6),
    axis.title.y = element_text(size= 6),
    legend.position = c(1.2,-0.1),
    legend.justification = c("right", "center"),
    legend.box.just = "right",
    legend.direction = "vertical"
  )

# Save the plot
ggsave(plot,filename=paste0("ADO_plot_summary.png"),width = 4,height = 3)

# Filter the 20-80% interval from the main df
df_filtered <- ado %>%
  filter(!grepl("0-0.1|0.1-0.2|0.8-0.9|0.9-1", File_Interval))

# Samplewise summary of ADO_20-80%
df_summary20_80 <- df_filtered %>%
  mutate(sample_name = sub("_[^_]+$", "", File_Interval)) %>%
  group_by(sample_name) %>%
  summarize(ADO_PERC = sum(Prop))

# Save the summary table
write.table(df_summary20_80,file = paste0("merged_ADO_summary.tsv"),
        append = F,sep = "\t",col.names = T,row.names = F, quote =F)
