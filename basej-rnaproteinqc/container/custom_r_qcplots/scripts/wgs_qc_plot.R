library(ohchibi)
library(ggpubr)
library(dplyr)
library(optparse)

# Shared 1-5 tier scoring, and the per-assay threshold specs.
# Installed to /usr/local/bin in the container; falls back to the script's own directory so
# the script can be run straight from the repo.
for (.f in c("qc_scoring.R", "qc_scoring_specs.R")) {
  .p <- c(file.path("/usr/local/bin", .f), .f)
  .hit <- .p[file.exists(.p)]
  if (!length(.hit)) stop(sprintf("%s not found in /usr/local/bin or the working directory", .f))
  source(.hit[1])
}

set.seed(130816)

bskb_col<-c("#12284C", "#1082A2","#A0CC2C" , "#DD14D3","#F45D34","#777776", "#FFFFFF")


# Established WGS gates; tier 4 of the scoring uses exactly these values, so tier >= 4 is the
# same verdict as the previous "all five pass".
cutoff_num_reads <- 50000000
cutoff_pct_dup<- 0.25
cutoff_pct_chim<- 0.15
cutoff_1x <- 0.9
cutoff_5x <- 0.7

# Tier 5 (Excellent): tighter coverage breadth and duplication. Chimeras stay at the gate
# value because on real data the chimera gate is almost never the binding constraint.
tight_pct_dup <- 0.15
tight_1x <- 0.95
tight_5x <- 0.85

# Tier 3 (Borderline): usable if capacity allows.
borderline_reads <- 25000000
borderline_1x <- 0.8
borderline_5x <- 0.5

# Below this many reads the coverage metrics measure depth rather than the library, so the
# sample is reported Inconclusive instead of failed. Expressed as a fraction of the read
# target so it tracks any retuning of cutoff_num_reads.
depth_floor_frac_of_target <- 0.1


plot_qc_wgs <- function(metrics_file) {

df <- read.table(file = metrics_file,header = TRUE,sep = "\t") %>%
  dplyr::rename(.data =.,SampleId = biosample)



Tab <- df[,c("total_reads","pct_duplication","pct_chimeras","pct_1x","pct_5x")]
rownames(Tab) <- df$SampleId

#Cluster samples pattern
mclust_samples <- hclust(d = dist(Tab),method = "ward.D")

order_samples <- mclust_samples$order %>% mclust_samples$labels[.]

df$SampleId <- df$SampleId %>% factor(levels = order_samples)

#### Assign the 1-5 quality tier (see qc_scoring.R)
#
# Replaces the previous count of passed thresholds. A tier requires ALL of its conditions, so
# a sample cannot rank above its weakest metric - under the count, missing the read target by
# 5% cost the same single point as duplication being three times over, and a sample failing a
# gate outright still rendered in the second band. Tier 4 keeps the gate values, so tier >= 4
# is identical to the old "all five pass".
df_tiers <- qc_assign_tiers(
  df,
  tiers = qc_spec_wgs(cutoff_num_reads = cutoff_num_reads,
                      cutoff_pct_dup = cutoff_pct_dup,
                      cutoff_pct_chim = cutoff_pct_chim,
                      cutoff_1x = cutoff_1x,
                      cutoff_5x = cutoff_5x,
                      tight_pct_dup = tight_pct_dup,
                      tight_1x = tight_1x,
                      tight_5x = tight_5x,
                      borderline_reads = borderline_reads,
                      borderline_1x = borderline_1x,
                      borderline_5x = borderline_5x),
  below_floor = qc_depth_floor(df$total_reads,
                               min_value = depth_floor_frac_of_target * cutoff_num_reads)
)

df_clust <- data.frame(SampleId = df_tiers$SampleId,
                       Cluster = df_tiers$TierLabel,
                       Band = df_tiers$Band,
                       stringsAsFactors = FALSE)

#Append tier information

df$Cluster <- match(df$SampleId,df_clust$SampleId) %>% df_clust$Cluster[.]

#### Total reads ####
df_ag <- df[,c("SampleId","Cluster","total_reads")] %>% unique %>%
  aggregate(total_reads~Cluster,.,quantile)

colnames(df_ag$total_reads) <- colnames(df_ag$total_reads) %>%
  gsub(pattern = "%",replacement = "") %>%
  gsub(pattern = "^",replacement = "PercTotReads")

df <- match(df$Cluster,df_ag$Cluster) %>% df_ag[.,-1] %>% cbind(df,.)


p_totreads <- ggplot(data = df,aes(total_reads,SampleId)) +
  geom_rect(aes(xmin =PercTotReads25 ,xmax =PercTotReads75,ymin = -Inf,ymax = Inf,group = Cluster),
            fill = "#A0CC2C",color = NA,alpha = 0.1)+
  geom_vline(aes(xintercept = PercTotReads50, group = Cluster), colour = '#DD14D3',
             size = 1.3) +
  geom_vline(xintercept = cutoff_num_reads,color = "red",linetype = "longdash")+
  geom_line(group = 1)+
  #geom_bar(stat = "identity",,width = 1,fill = "black",color = NA) +
  facet_grid(Cluster~.,space = "free",scales = "free") +
  theme_ohchibi(size_panel_border = 0.3)+
  theme(
    legend.position = "top",
    panel.grid.major.x = element_line(linetype = "dotted",color= "grey"),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_line(size = unit(0.1,"line")),
    axis.title.y = element_blank(),
    axis.ticks.x = element_line(size = unit(0.1,"line")),
    axis.text.x = element_text(size = 9,color = "grey30",angle = 90,vjust = 0.5,hjust = 1),
    panel.background = element_blank(),
    #panel.border = element_blank(),
    axis.line = element_line(color = 'black',size = 0.3),
    strip.text.y = element_blank(),
    panel.spacing.y = unit(0.1, "lines"),
    axis.title.x = element_text(size = 9)
  ) +
  xlab(label = "Total number of input reads")  +
  scale_x_log10()


#### PCT DUP####
df_ag <- df[,c("SampleId","Cluster","pct_duplication")] %>% unique %>%
  aggregate(pct_duplication~Cluster,.,quantile)

colnames(df_ag$pct_duplication) <- colnames(df_ag$pct_duplication) %>%
  gsub(pattern = "%",replacement = "") %>%
  gsub(pattern = "^",replacement = "PercPCT_DUPLICATION")

df <- match(df$Cluster,df_ag$Cluster) %>% df_ag[.,-1] %>% cbind(df,.)


p_dup <- ggplot(data = df,aes(pct_duplication,SampleId)) +
  geom_rect(aes(xmin =PercPCT_DUPLICATION25 ,xmax =PercPCT_DUPLICATION75,ymin = -Inf,ymax = Inf,group = Cluster),
            fill = "#A0CC2C",color = NA,alpha = 0.1)+
  geom_vline(aes(xintercept = PercPCT_DUPLICATION50, group = Cluster), colour = '#DD14D3',
             size = 1.3) +
  geom_vline(xintercept = cutoff_pct_dup,color = "red",linetype = "longdash")+
  geom_line(group = 1)+
  #geom_bar(stat = "identity",,width = 1,fill = "black",color = NA) +
  facet_grid(Cluster~.,space = "free",scales = "free") +
  theme_ohchibi(size_panel_border = 0.3)+
  theme(
    legend.position = "top",
    panel.grid.major.x = element_line(linetype = "dotted",color= "grey"),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.x = element_line(size = unit(0.1,"line")),
    axis.text.x = element_text(size = 9,color = "grey30",angle = 90,vjust = 0.5,hjust = 1),
    panel.background = element_blank(),
    #panel.border = element_blank(),
    axis.line = element_line(color = 'black',size = 0.3),
    strip.text.y = element_blank(),
    panel.spacing.y = unit(0.1, "lines"),
    axis.title.x = element_text(size = 9)
  ) +
  xlab(label = "pct_duplication")  +
  scale_x_continuous(breaks = seq(0,1,0.1),limits = c(0,1))


#### PCT CHIM####
df_ag <- df[,c("SampleId","Cluster","pct_chimeras")] %>% unique %>%
  aggregate(pct_chimeras~Cluster,.,quantile)

colnames(df_ag$pct_chimeras) <- colnames(df_ag$pct_chimeras) %>%
  gsub(pattern = "%",replacement = "") %>%
  gsub(pattern = "^",replacement = "PercPCT_CHIMERAS")

df <- match(df$Cluster,df_ag$Cluster) %>% df_ag[.,-1] %>% cbind(df,.)


p_chim <- ggplot(data = df,aes(pct_chimeras,SampleId)) +
  geom_rect(aes(xmin =PercPCT_CHIMERAS25 ,xmax =PercPCT_CHIMERAS75,ymin = -Inf,ymax = Inf,group = Cluster),
            fill = "#A0CC2C",color = NA,alpha = 0.1)+
  geom_vline(aes(xintercept = PercPCT_CHIMERAS50, group = Cluster), colour = '#DD14D3',
             size = 1.3) +
  geom_vline(xintercept = cutoff_pct_chim,color = "red",linetype = "longdash")+
  geom_line(group = 1)+
  #geom_bar(stat = "identity",,width = 1,fill = "black",color = NA) +
  facet_grid(Cluster~.,space = "free",scales = "free") +
  theme_ohchibi(size_panel_border = 0.3)+
  theme(
    legend.position = "top",
    panel.grid.major.x = element_line(linetype = "dotted",color= "grey"),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.x = element_line(size = unit(0.1,"line")),
    axis.text.x = element_text(size = 9,color = "grey30",angle = 90,vjust = 0.5,hjust = 1),
    panel.background = element_blank(),
    #panel.border = element_blank(),
    axis.line = element_line(color = 'black',size = 0.3),
    strip.text.y = element_blank(),
    panel.spacing.y = unit(0.1, "lines"),
    axis.title.x = element_text(size = 9)
  ) +
  xlab(label = "pct_chimeras")  +
  scale_x_continuous(breaks = seq(0,1,0.1),limits = c(0,1))


#### PCT 1X####
df_ag <- df[,c("SampleId","Cluster","pct_1x")] %>% unique %>%
  aggregate(pct_1x~Cluster,.,quantile)

colnames(df_ag$pct_1x) <- colnames(df_ag$pct_1x) %>%
  gsub(pattern = "%",replacement = "") %>%
  gsub(pattern = "^",replacement = "PercPCT_1X")

df <- match(df$Cluster,df_ag$Cluster) %>% df_ag[.,-1] %>% cbind(df,.)


p_1x <- ggplot(data = df,aes(pct_1x,SampleId)) +
  geom_rect(aes(xmin =PercPCT_1X25 ,xmax =PercPCT_1X75,ymin = -Inf,ymax = Inf,group = Cluster),
            fill = "#A0CC2C",color = NA,alpha = 0.1)+
  geom_vline(aes(xintercept = PercPCT_1X50, group = Cluster), colour = '#DD14D3',
             size = 1.3) +
  geom_vline(xintercept = cutoff_1x,color = "red",linetype = "longdash")+
  geom_line(group = 1)+
  #geom_bar(stat = "identity",,width = 1,fill = "black",color = NA) +
  facet_grid(Cluster~.,space = "free",scales = "free") +
  theme_ohchibi(size_panel_border = 0.3)+
  theme(
    legend.position = "top",
    panel.grid.major.x = element_line(linetype = "dotted",color= "grey"),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.x = element_line(size = unit(0.1,"line")),
    axis.text.x = element_text(size = 9,color = "grey30",angle = 90,vjust = 0.5,hjust = 1),
    panel.background = element_blank(),
    #panel.border = element_blank(),
    axis.line = element_line(color = 'black',size = 0.3),
    strip.text.y = element_blank(),
    panel.spacing.y = unit(0.1, "lines"),
    axis.title.x = element_text(size = 9)
  ) +
  xlab(label = "pct_1x")  +
  scale_x_continuous(breaks = seq(0,1,0.1),limits = c(0,1))


#### PCT 5X####
df_ag <- df[,c("SampleId","Cluster","pct_5x")] %>% unique %>%
  aggregate(pct_5x~Cluster,.,quantile)

colnames(df_ag$pct_5x) <- colnames(df_ag$pct_5x) %>%
  gsub(pattern = "%",replacement = "") %>%
  gsub(pattern = "^",replacement = "PercPCT_5X")

df <- match(df$Cluster,df_ag$Cluster) %>% df_ag[.,-1] %>% cbind(df,.)


p_5x <- ggplot(data = df,aes(pct_5x,SampleId)) +
  geom_rect(aes(xmin =PercPCT_5X25 ,xmax =PercPCT_5X75,ymin = -Inf,ymax = Inf,group = Cluster),
            fill = "#A0CC2C",color = NA,alpha = 0.1)+
  geom_vline(aes(xintercept = PercPCT_5X50, group = Cluster), colour = '#DD14D3',
             size = 1.3) +
  geom_vline(xintercept = cutoff_5x,color = "red",linetype = "longdash")+
  geom_line(group = 1)+
  #geom_bar(stat = "identity",,width = 1,fill = "black",color = NA) +
  facet_grid(Cluster~.,space = "free",scales = "free") +
  theme_ohchibi(size_panel_border = 0.3)+
  theme(
    legend.position = "top",
    panel.grid.major.x = element_line(linetype = "dotted",color= "grey"),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.x = element_line(size = unit(0.1,"line")),
    axis.text.x = element_text(size = 9,color = "grey30",angle = 90,vjust = 0.5,hjust = 1),
    panel.background = element_blank(),
    #panel.border = element_blank(),
    axis.line = element_line(color = 'black',size = 0.3),
    strip.text.y = element_text(size = 9,angle = 0),
    panel.spacing.y = unit(0.1, "lines"),
    axis.title.x = element_text(size = 9)
  ) +
  xlab(label = "pct_5x")  +
  scale_x_continuous(breaks = seq(0,1,0.1),limits = c(0,1))




#Try another approach to visualize the results

melted <- colnames(df) %>% grep(pattern = "pct_.*x$",value = F) %>% c(1,.) %>% df[,.] %>%
  reshape2::melt() %>%
  dplyr::mutate(.data =.,variable = variable %>% gsub(pattern = "pct_",replacement = ""))


options("scipen"=100, "digits"=10)
levels_fac <- paste0(melted$variable %>% gsub(pattern = "x",replacement = "")  %>% unique %>% as.numeric %>% sort,"x")


melted$variable <- melted$variable %>% factor(levels = levels_fac)
melted$Cluster <- match(melted$SampleId,df$SampleId) %>% df$Cluster[.]


paleta_cluster <- bskb_col[1:6] %>% rev
names(paleta_cluster) <- c("5","4","3","2","1","0")


p1 <- ggplot(data =melted,aes(variable,value)) +
  geom_line(aes(group = SampleId,color = Cluster))+
  geom_point()  +
  stat_summary(fun = median,geom = "point",color = "red",shape = 15,size = 2)+
  #facet_grid(.~Cluster,space = "free",scales = "free") +
  scale_y_continuous(breaks = seq(0,1,0.1),limits = c(0,1)) +
  theme_ohchibi(size_panel_border = 0.3,size_title_text = 10,size_legend_text = 9,size_axis_title.x = 10,size_axis_title.y = 10) +
  theme(
    legend.position = "top",
    panel.grid.major.x = element_line(linetype = "dotted",color = "grey"),
    panel.grid.major.y = element_line(linetype = "dotted",color = "grey"),
    axis.text.y = element_text(size = 10,color = "grey30"),
    axis.title.y = element_text(size = 10,color = "black"),
    axis.ticks.y = element_line(size = unit(0.1,"line")),
    axis.ticks.x = element_line(size = unit(0.1,"line")),
    axis.text.x = element_text(size = 10,color = "grey30",angle = 90,vjust = 0.5,hjust = 1),
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = 'black',size = 0.3),
    strip.text.y = element_blank(),
    panel.spacing.y = unit(0, "lines"),
    axis.title.x = element_text(size = 10,color = "black")) +
  xlab(label = "Coverage") +
  ylab(label = "Proportion of total positions") +
  scale_color_manual(values = paleta_cluster)


grDevices::pdf(file = NULL)
composition1 <- egg::ggarrange(p_totreads,p_dup,p_chim,p_1x,p_5x,nrow = 1)
grDevices::dev.off()

#Send table with verdict
df <- df %>%
  dplyr::relocate(.data =.,c("SampleId","Cluster"))  %>%
  dplyr::arrange(.data =.,Cluster)
p2 <- df$Cluster %>% table %>%
  data.frame %>%
  dplyr::rename(.data =.,QualityCluster = ".",NumCells = Freq) %>%
  dplyr::mutate(.data =.,PropCells =round( NumCells/nrow(df),3)) %>%
  gridExtra::tableGrob(rows = NULL) %>%
  as_ggplot() +
  theme(
    plot.title = element_text(size = 15,vjust =0,hjust = 0.5)
  )+
  ggtitle(label = "Total of usable cells")

grDevices::pdf(file = NULL)
composition2 <- egg::ggarrange(p1,p2,nrow = 1,widths = c(1,0.6),labels = c("B","C"))
grDevices::dev.off()

composition_qc_wgs <- cowplot::plot_grid(composition1,composition2,nrow = 2,labels = c(" "," "),label_y = 1,rel_heights = c(1,0.75))


title <- cowplot::ggdraw() +
  cowplot::draw_label("Quality control of WGS dataset",
                      size = 13, x = 0.5, vjust = 0)


final_plot <- cowplot::plot_grid(title, composition_qc_wgs, ncol = 1, rel_heights = c(0.05, 1))



# `CompositeScore` is retained as the numeric 1-5 tier so downstream consumers and the
# published WGS-QC_ConsensusScores.txt keep a stable column name; QC_Band is the
# Pass/Borderline/Fail collapse and BlockingMetric says which condition held each sample back.
merged <- df_tiers %>%
  dplyr::transmute(SampleId = as.character(SampleId),
                   CompositeScore = Tier,
                   QC_Label = as.character(TierLabel),
                   QC_Band = as.character(Band),
                   BlockingMetric = BlockingMetric) %>%
  merge(df %>% dplyr::mutate(.data =., SampleId = as.character(SampleId)) %>%
          dplyr::select(dplyr::any_of(
            c("SampleId","total_reads","pct_duplication","pct_chimeras","pct_1x","pct_5x"))) %>%
          unique,
        by = "SampleId", all.x = TRUE) %>%
  dplyr::arrange(.data =., dplyr::desc(CompositeScore), SampleId)

df_sum_cat <- merged$QC_Label %>% factor(levels = QC_TIER_ORDER) %>%
  table %>% data.frame %>%
  dplyr::rename(.data =., Category = ".", NumberSamples = Freq) %>%
  dplyr::mutate(.data =., ProportionSamples = ((NumberSamples/sum(NumberSamples))*100) %>% round(2))

df_sum_band <- merged$QC_Band %>% factor(levels = QC_BAND_ORDER) %>%
  table %>% data.frame %>%
  dplyr::rename(.data =., Band = ".", NumberSamples = Freq) %>%
  dplyr::mutate(.data =., ProportionSamples = ((NumberSamples/sum(NumberSamples))*100) %>% round(2))

return(list(df_verdict = merged,
            df_sum_verdict = df_sum_cat,
            df_sum_band = df_sum_band,
            composition = final_plot))
}


option_list <- list(
  make_option("--metrics_file", type = "character",
              help = "WGS metrics TSV file (must contain biosample column)")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$metrics_file)) {
  stop("Please provide --metrics_file")
}

res <- plot_qc_wgs(metrics_file = opt$metrics_file)

ggplot2::ggsave(filename = "qc_wgs.pdf", plot = res$composition,
                width = 12, height = 14, units = "in")

ggplot2::ggsave(filename = "WGS-QC_composition_mqc.jpg", plot = res$composition,
                width = 12, height = 14, units = "in", dpi = 300)

write.table(res$df_verdict, "WGS-QC_ConsensusScores.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(res$df_sum_verdict, "WGS-QC_ConsensusScores_SummaryTable_mqc.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(res$df_sum_band, "WGS-QC_QCBand_SummaryTable_mqc.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
