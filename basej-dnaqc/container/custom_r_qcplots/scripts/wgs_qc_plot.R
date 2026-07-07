library(ohchibi)
library(ggpubr)
library(dplyr)
library(optparse)


set.seed(130816)

bskb_col<-c("#12284C", "#1082A2","#A0CC2C" , "#DD14D3","#F45D34","#777776", "#FFFFFF")


cutoff_num_reads <- 50000000
cutoff_pct_dup<- 0.25
cutoff_pct_chim<- 0.15
cutoff_1x <- 0.9
cutoff_5x <- 0.7


plot_qc_wgs <- function(metrics_file) {

df <- read.table(file = metrics_file,header = TRUE,sep = "\t") %>%
  dplyr::rename(.data =.,SampleId = biosample)



Tab <- df[,c("total_reads","pct_duplication","pct_chimeras","pct_1x","pct_5x")]
rownames(Tab) <- df$SampleId

#Cluster samples pattern
mclust_samples <- hclust(d = dist(Tab),method = "ward.D")

order_samples <- mclust_samples$order %>% mclust_samples$labels[.]

df$SampleId <- df$SampleId %>% factor(levels = order_samples)

#Quantify verdict samples
#### Construct the clustered version of the figures with composite score logic
Mat <- matrix(data = 0,nrow = nrow(df),ncol = 5)
rownames(Mat) <- df$SampleId
Mat[,1][which(df$total_reads >= cutoff_num_reads)] <- 1
Mat[,2][which(df$pct_duplication < cutoff_pct_dup)] <- 1
Mat[,3][which(df$pct_chimeras < cutoff_pct_chim)] <- 1
Mat[,4][which(df$pct_1x >= cutoff_1x)] <- 1
Mat[,5][which(df$pct_5x >= cutoff_5x)] <- 1


#Create sum and that will become the cluster
df_clust <- rowSums(Mat) %>% data.frame %>%
  tibble::rownames_to_column() %>%
  dplyr::rename(.data =.,SampleId = rowname, Cluster = ".")

df_clust$Cluster <- df_clust$Cluster %>% factor(levels = c("5","4","3","2","1","0"))

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



merged <- Mat %>% as.data.frame %>%
  tibble::rownames_to_column() %>%
  dplyr::rename(.data =., SampleId = rowname,
                Verdict_TOTAL_READS = V1,
                Verdict_PCT_DUPLICATION = V2,
                Verdict_PCT_CHIMERAS = V3,
                Verdict_PCT_1X = V4,
                Verdict_PCT_5X = V5) %>%
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
