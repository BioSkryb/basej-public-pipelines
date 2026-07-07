library(ohchibi)
library(dplyr)
library(optparse)

set.seed(130816)

bskb_col<-c("#12284C", "#1082A2","#A0CC2C" , "#DD14D3","#F45D34","#777776", "#FFFFFF")


cutoff_num_reads <- 500000
cutoff_10x <- 0.75
cutoff_zero_cov <- 0.05
cutoff_fold_80 <- 5


plot_qc_wes <- function(metrics_file) {

df <- read.table(file = metrics_file,header = TRUE,sep = "\t") %>%
  dplyr::rename(.data =.,SampleId = biosample)

df$fold_80_base_penalty <- df$fold_80_base_penalty %>% as.numeric


Tab <- df[,c("total_reads","pct_target_bases_10x","zero_cvg_targets_pct","fold_80_base_penalty")]
rownames(Tab) <- df$SampleId

#Cluster samples pattern
mclust_samples <- hclust(d = dist(Tab),method = "ward.D")

order_samples <- mclust_samples$order %>% mclust_samples$labels[.]

df$SampleId <- df$SampleId %>% factor(levels = order_samples)

#Quantify verdict samples
#### Construct the clustered version of the figures with composite score logic
Mat <- matrix(data = 0,nrow = nrow(df),ncol = 4)
rownames(Mat) <- df$SampleId
Mat[,1][which(df$total_reads >= cutoff_num_reads)] <- 1
Mat[,2][which(df$pct_target_bases_10x >= cutoff_10x)] <- 1
Mat[,3][which(df$zero_cvg_targets_pct < cutoff_zero_cov)] <- 1
Mat[,4][which(df$fold_80_base_penalty < cutoff_fold_80)] <- 1

#Create sum and that will become the cluster
df_clust <- rowSums(Mat) %>% data.frame %>%
  tibble::rownames_to_column() %>%
  dplyr::rename(.data =.,SampleId = rowname, Cluster = ".")

df_clust$Cluster <- df_clust$Cluster %>% factor(levels = c("4","3","2","1","0"))

#Append cluster information

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


#### PCT 10 x ####
df_ag <- df[,c("SampleId","Cluster","pct_target_bases_10x")] %>% unique %>%
  aggregate(pct_target_bases_10x~Cluster,.,quantile)

colnames(df_ag$pct_target_bases_10x) <- colnames(df_ag$pct_target_bases_10x) %>%
  gsub(pattern = "%",replacement = "") %>%
  gsub(pattern = "^",replacement = "PercTARGET_BASES_10X")

df <- match(df$Cluster,df_ag$Cluster) %>% df_ag[.,-1] %>% cbind(df,.)


p_10x <- ggplot(data = df,aes(pct_target_bases_10x,SampleId)) +
  geom_rect(aes(xmin =PercTARGET_BASES_10X25 ,xmax =PercTARGET_BASES_10X75,ymin = -Inf,ymax = Inf,group = Cluster),
            fill = "#A0CC2C",color = NA,alpha = 0.1)+
  geom_vline(aes(xintercept = PercTARGET_BASES_10X50, group = Cluster), colour = '#DD14D3',
             size = 1.3) +
  geom_vline(xintercept = cutoff_10x,color = "red",linetype = "longdash")+
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
  xlab(label = "pct_target_bases_10x")  +
  scale_x_continuous(breaks = seq(0,1,0.1),limits = c(0,1))



#### Zero coverage ####
df_ag <- df[,c("SampleId","Cluster","zero_cvg_targets_pct")] %>% unique %>%
  aggregate(zero_cvg_targets_pct~Cluster,.,quantile)

colnames(df_ag$zero_cvg_targets_pct) <- colnames(df_ag$zero_cvg_targets_pct) %>%
  gsub(pattern = "%",replacement = "") %>%
  gsub(pattern = "^",replacement = "PercZERO_CVG_TARGETS_PCT")

df <- match(df$Cluster,df_ag$Cluster) %>% df_ag[.,-1] %>% cbind(df,.)


p_zero <- ggplot(data = df,aes(zero_cvg_targets_pct,SampleId)) +
  geom_rect(aes(xmin =PercZERO_CVG_TARGETS_PCT25 ,xmax =PercZERO_CVG_TARGETS_PCT75,ymin = -Inf,ymax = Inf,group = Cluster),
            fill = "#A0CC2C",color = NA,alpha = 0.1)+
  geom_vline(aes(xintercept = PercZERO_CVG_TARGETS_PCT50, group = Cluster), colour = '#DD14D3',
             size = 1.3) +
  geom_vline(xintercept = cutoff_zero_cov,color = "red",linetype = "longdash")+
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
  xlab(label = "zero_cvg_targets_pct")  +
  scale_x_continuous(breaks = seq(0,1,0.1),limits = c(0,1))



####  Fold 80 base ####
df_ag <- df[,c("SampleId","Cluster","fold_80_base_penalty")] %>% unique %>%
  aggregate(fold_80_base_penalty~Cluster,.,quantile)

colnames(df_ag$fold_80_base_penalty) <- colnames(df_ag$fold_80_base_penalty) %>%
  gsub(pattern = "%",replacement = "") %>%
  gsub(pattern = "^",replacement = "PercFOLD_80_BASE_PENALTY")

df <- match(df$Cluster,df_ag$Cluster) %>% df_ag[.,-1] %>% cbind(df,.)


p_80 <- ggplot(data = df,aes(fold_80_base_penalty,SampleId)) +
  geom_rect(aes(xmin =PercFOLD_80_BASE_PENALTY25 ,xmax =PercFOLD_80_BASE_PENALTY75,ymin = -Inf,ymax = Inf,group = Cluster),
            fill = "#A0CC2C",color = NA,alpha = 0.1)+
  geom_vline(aes(xintercept = PercFOLD_80_BASE_PENALTY50, group = Cluster), colour = '#DD14D3',
             size = 1.3) +
  geom_vline(xintercept = cutoff_fold_80,color = "red",linetype = "longdash")+
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
  xlab(label = "fold_80_base_penalty") +
  scale_x_continuous(breaks = seq(0,10,1),limits = c(0,10),oob = squish)



#Try another approach to visualize the results

melted <- colnames(df) %>% grep(pattern = "pct_target_bases_",value = F) %>% c(1,.) %>% df[,.] %>%
  reshape2::melt() %>%
  dplyr::mutate(.data =.,variable = variable %>% gsub(pattern = "pct_target_bases_",replacement = ""))


options("scipen"=100, "digits"=10)
levels_fac <- paste0(melted$variable %>% gsub(pattern = "x",replacement = "")  %>% unique %>% as.numeric %>% sort,"x")


melted$variable <- melted$variable %>% factor(levels = levels_fac)
melted$Cluster <- match(melted$SampleId,df$SampleId) %>% df$Cluster[.]


paleta_cluster <- bskb_col[1:5] %>% rev
names(paleta_cluster) <- c("4","3","2","1","0")


p1 <- ggplot(data =melted,aes(variable,value)) +
  geom_line(aes(group = SampleId,color = Cluster))+
  geom_point()  +
  stat_summary(fun = median,geom = "point",color = "red",shape = 15,size = 2)+
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

#Plot about fold coverage

df$PROP_GOOD_BASES <- (df$on_bait_bases+df$near_bait_bases)/(df$on_bait_bases+df$near_bait_bases+df$off_bait_bases)
df$PROP_BAD_BASES<- (df$off_bait_bases)/(df$on_bait_bases+df$near_bait_bases+df$off_bait_bases)

#Now load phylogeny and label using the metrics output
melted <- reshape2::melt(data = df,id.vars = c("SampleId"),measure.vars = c("PROP_GOOD_BASES","PROP_BAD_BASES"))


melted$variable <- melted$variable %>%
  gsub(pattern = "PROP_GOOD_BASES",replacement= "ON_BAIT_BASES+NEAR_BAIT_BASES") %>%
  gsub(pattern = "PROP_BAD_BASES",replacement= "OFF_BAIT_BASES")  %>%
  factor(levels = c("ON_BAIT_BASES+NEAR_BAIT_BASES","OFF_BAIT_BASES") %>% rev)


melted$Cluster <- match(melted$SampleId,df$SampleId) %>% df$Cluster[.]

p3 <- ggplot(data = melted %>% subset(Cluster != 0),aes(SampleId,value)) +
  geom_bar(stat = "identity",aes(fill = variable),color = "transparent",width = 1) +
  facet_grid(.~Cluster,space = "free",scales = "free") +
  theme_ohchibi(size_panel_border = 0.3,size_title_text = 10,size_legend_text = 9,size_axis_title.x = 10,size_axis_title.y = 10) +
  theme(
    legend.position = "top",
    panel.grid.major.x = element_line(linetype = "dotted",color = "grey"),
    panel.grid.major.y = element_line(linetype = "dotted",color = "grey"),
    axis.text.y = element_text(size = 10,color = "grey30"),
    axis.title.y = element_text(size = 10,color = "black"),
    axis.ticks.y = element_line(size = unit(0.1,"line")),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = 'black',size = 0.3),
    strip.text.y = element_blank(),
    panel.spacing.y = unit(0, "lines"),
    axis.title.x = element_text(size = 10,color = "black")) +
  scale_y_continuous(breaks = seq(0,1,0.1),expand = c(0,0)) +
  scale_fill_manual(values = c("#12284C","#A0CC2C")  %>% rev,name = "Category") +
  xlab(label = "Single Cells") +
  ylab(label = "Proportion of total bases")

#Check the other metrics

grDevices::pdf(file = NULL)
composition1 <- egg::ggarrange(p_totreads,p_10x,p_zero,p_80,nrow = 1)
grDevices::dev.off()



composition2 <- ggpubr::ggarrange(p1,p3,nrow = 1,labels = c("B","C"))

composition_qc_wes <- cowplot::plot_grid(composition1,composition2,nrow = 2,labels = c("",""),label_y = 1,rel_heights = c(1,0.75))


title <- cowplot::ggdraw() +
  cowplot::draw_label("Quality control of WES dataset",
                      size = 16, x = 0.5, vjust = 0)

final_plot <- cowplot::plot_grid(title, composition_qc_wes, ncol = 1, rel_heights = c(0.05, 1))

merged <- Mat %>% as.data.frame %>%
  tibble::rownames_to_column() %>%
  dplyr::rename(.data =., SampleId = rowname,
                Verdict_TOTAL_READS = V1,
                Verdict_PCT_TARGET_BASES_10X = V2,
                Verdict_ZERO_CVG_TARGETS_PCT = V3,
                Verdict_FOLD_80_BASE_PENALTY = V4) %>%
  merge(df_clust, ., by = "SampleId") %>%
  dplyr::rename(.data =., CompositeScore = Cluster)

df_sum_cat <- merged$CompositeScore %>% table %>% data.frame %>%
  dplyr::rename(.data =., Category = ".", NumberSamples = Freq) %>%
  dplyr::mutate(.data =., ProportionSamples = ((NumberSamples/sum(NumberSamples))*100) %>% round(2))

return(list(df_verdict = merged,
            df_sum_verdict = df_sum_cat,
            composition = final_plot))
}


option_list <- list(
  make_option("--metrics_file", type = "character",
              help = "WES metrics TSV file (must contain biosample column)")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$metrics_file)) {
  stop("Please provide --metrics_file")
}

res <- plot_qc_wes(metrics_file = opt$metrics_file)

ggplot2::ggsave(filename = "qc_wes.pdf", plot = res$composition,
                width = 12, height = 14, units = "in")

ggplot2::ggsave(filename = "WES-QC_composition_mqc.jpg", plot = res$composition,
                width = 12, height = 14, units = "in", dpi = 300)

write.table(res$df_verdict, "WES-QC_ConsensusScores.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(res$df_sum_verdict, "WES-QC_ConsensusScores_SummaryTable_mqc.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
