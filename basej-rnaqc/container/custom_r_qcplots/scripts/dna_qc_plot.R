library(ohchibi)
library(dplyr)
library(ggtree)
library(optparse)
library(jsonlite)

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

plot_qc_dna <- function(
    seg_copy_file= NULL,
    metrics_file =  NULL,
    metadata_file = NULL,
    cutoff_preseq = 3500000000,
    cutoff_chimeras = 0.2,
    cutoff_chrM  = 0.2,
    cutoff_cnv_mapd = 0.25,
    cutoff_cnv_sk = 0.25,
    # tier-5 (Excellent) tightening of the CNV metrics; tier 4 keeps the gate values above
    tight_cnv_mapd = 0.20,
    tight_cnv_sk = 0.20,
    # below this many reads the coverage metrics measure depth rather than the cell,
    # so the cell is reported Inconclusive instead of failed
    cutoff_reads_floor = 1e6,
    size_labels_samples = 1,
    relative_size_cnv_plot = 1,
    relative_size_preseq_plot = 0.2,
    relative_size_mapd_plot = 0.2,
    relative_size_skew_plot = 0.2,
    relative_size_chrm_plot = 0.2,
    relative_size_chimeras_plot = 0.2,
    relative_size_bar_plot_group = 0.03
    
){
  
  #Read tables and determine usable samples
  df_cnv <- read.table(file = seg_copy_file,header = TRUE,check.names = FALSE)
  ids_cnv <- colnames(df_cnv)[4:ncol(df_cnv)]
  
  df_preq <- read.table(file = metrics_file,header = TRUE,sep = "\t") %>%
    dplyr::rename(.data =.,SampleId = sample_name)
  ids_metrics <- df_preq$SampleId
  
  #Check if metadata file was passed
  if(is.null(metadata_file)){
    
    usable_ids <- intersect(ids_cnv,ids_metrics)
    df_meta <- NULL
    
  }else{
    
    df_meta <- read.table(file = metadata_file,header = TRUE,sep = ",")
    
    # Handle group/groups column FIRST before any other operations
    # This prevents duplicate column name issues downstream
    if ("group" %in% colnames(df_meta) && "groups" %in% colnames(df_meta)) {
      # Both exist - remove 'groups' and keep 'group'
      df_meta <- df_meta[, colnames(df_meta) != "groups", drop = FALSE]
    } else if ("groups" %in% colnames(df_meta) && !("group" %in% colnames(df_meta))) {
      # Only 'groups' exists - rename it to 'group'
      colnames(df_meta)[colnames(df_meta) == "groups"] <- "group"
    } else if (!("group" %in% colnames(df_meta)) && !("groups" %in% colnames(df_meta))) {
      # Neither exists - create default group for all samples
      df_meta$group <- "Group1"
    }
    # If only 'group' exists, no action needed
    
    # Now handle biosampleName column
    colnames(df_meta)[colnames(df_meta) == "biosampleName"] <- "SampleId"
    
    # Filter to only samples with non-null and non-blank group information
    df_meta <- df_meta[!is.na(df_meta$group) & df_meta$group != "", ]
    ids_meta <- df_meta$SampleId
    usable_ids <- intersect(ids_cnv,ids_metrics) %>% intersect(ids_meta)
    
  }
  
  #### Processs CNV results ####
  melted <- reshape2::melt(data = df_cnv,id.vars = 1:3) %>%
    dplyr::rename(.data = .,SampleId = variable) %>%
    dplyr::mutate(.data = .,SampleId = SampleId %>% gsub(pattern = "_sorted",replacement = ""))
  
  #Filter and keep only usable ids
  melted <- melted %>% 
    dplyr::filter(.data =.,SampleId %in% usable_ids) %>% droplevels
  
  melted$Range <- paste0(melted$START,"-",melted$END)
  
  #############################################################################
  ############# Modification to keep proper order of bins ####################
  melted$Range <- paste0(melted$START,"-",melted$END)
  order_range <- with(melted,order(START)) %>%
    melted[.,] %$% Range %>% unique
  melted$Range <- melted$Range %>% factor(levels = order_range)
  
  #melted$Range <- melted$Range %>% factor
  #############################################################################
  ##############################################################################
  
  # Dynamically determine chromosome order based on data present
  # Supports both human (chr1-22,X,Y) and mouse (chr1-19,X,Y) genomes
  all_chrs <- melted$CHR %>% unique %>% as.character
  # Extract numeric chromosomes and sort them numerically
  numeric_chrs <- all_chrs[grepl("^chr[0-9]+$", all_chrs)]
  numeric_order <- numeric_chrs[order(as.numeric(gsub("chr", "", numeric_chrs)))]
  # Add sex chromosomes in standard order if present
  sex_chrs <- c("chrX", "chrY")
  sex_chrs_present <- sex_chrs[sex_chrs %in% all_chrs]
  chr_levels <- c(numeric_order, sex_chrs_present)
  melted$CHR <- melted$CHR %>% factor(levels = chr_levels)


  #Determine blocks that do not change
  Tab <- acast(data = melted,formula = CHR+Range~SampleId,value.var = "value")
  # Ensure Tab is always a matrix (acast drops dimensions with single sample)
  if(is.null(dim(Tab))) {
    Tab <- matrix(Tab, ncol = 1, dimnames = list(names(Tab), unique(melted$SampleId)))
  }

  #Create a discrete palette
  paleta <- c("#4575b4","#abd9e9","#D9D9D9","#ffffbf","#fee090","#fdae61","#f46d43","#d73027","#a50026")
  melted$valueFactor <- melted$value
  melted$valueFactor[melted$valueFactor>=4] <- ">=4"
  melted$valueFactor <- melted$valueFactor %>% factor(levels = c("0","1","2","3",">=4"))

  paleta <- c("#313695","#74add1","#FAF9F6","#f46d43","#a50026")
  names(paleta) <- levels(melted$PloidyFactor)

  # Create chr column (without "chr" prefix) with dynamic ordering
  melted$chr <- melted$CHR %>%
    as.character %>%
    gsub(pattern = "chr",replacement = "")
  # Dynamically determine order based on data present
  all_chr_short <- melted$chr %>% unique
  numeric_chr_short <- all_chr_short[grepl("^[0-9]+$", all_chr_short)]
  numeric_order_short <- numeric_chr_short[order(as.numeric(numeric_chr_short))]
  sex_chr_short <- c("X", "Y")
  sex_chr_present <- sex_chr_short[sex_chr_short %in% all_chr_short]
  chr_levels_short <- c(numeric_order_short, sex_chr_present)
  melted$chr <- melted$chr %>% factor(levels = chr_levels_short)
  Tab_bray <- Tab
  Tab_bray[which(Tab_bray>8)] <- 8
  Tab_bray <- Tab_bray[which(rowSums(Tab_bray)  != 0), , drop = FALSE]
  
  # Handle single-sample case: hclust requires n >= 2

  if(ncol(Tab_bray) >= 2) {
    mdist <- vegdist(x = Tab_bray %>% t,method = "bray")
    mclust <- hclust(d = mdist,method = "ward.D")
    order_ids <- mclust$order %>% mclust$labels[.]
  } else {
    order_ids <- colnames(Tab_bray)
  }
  
  #Create ggtrew
  melted$SampleId <- melted$SampleId %>% factor(levels = order_ids)
  
  p2 <- ggplot(data = melted,aes(Range,SampleId)) +
    geom_raster(aes(fill = valueFactor)) +
    facet_grid(.~chr,space = "free",scales = "free") +
    theme_ohchibi(size_panel_border = 0.3) +
    scale_fill_manual(values = paleta ,name = "Estimated Ploidy") +
    theme(
      legend.position = "bottom",
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 1),
      axis.ticks.y = element_line(size = 0.2),
      strip.text.x = element_text(size = 9,color = "grey30"),
      strip.text.y = element_text(size = 9),
      legend.text = element_text(size = 9,color = "grey30"),
      legend.title = element_text(size = 10),
      plot.title = element_text(size = 10)
      
    ) +
    xlab(label = "Bin")  +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0))
  
  #Append dendrogram
  #p_tree_rows <- ggtree(mclust, ladderize = F, size = 0.3) +
  #  scale_y_continuous(expand =  c(0.001, 0.001))
  p_blank <- ggplot() + theme_void()
  
  
  #composition <- egg::ggarrange(p_tree_rows,p2,nrow =1 ,widths = c(0.1,1))
  
  #Check preseq values and plot them along
  df_preq <- df_preq %>% 
    dplyr::filter(.data =.,SampleId %in% usable_ids) %>% droplevels
  
  #Order them according to the CNV
  df_preq$SampleId <- df_preq$SampleId %>% factor(levels = order_ids)
  
  #### Assign the 1-5 quality tier (see qc_scoring.R)
  #
  # Replaces the previous count of passed thresholds. A tier now requires ALL of its
  # conditions, so a cell cannot rank above its weakest metric - under the old count, 53 of
  # 84 cells scoring 4 of 5 were actually failing PreSeq, 7 of them with PreSeq of zero.
  # Tier 4 keeps the original gate values, so tier >= 4 equals the old "all five pass".
  qc_require_metric(df_preq, c("MAPD_CNV_Log2", "SKEW_CNV"), what = "CNV")

  df_tiers <- qc_assign_tiers(
    df_preq,
    tiers = qc_spec_dna(cutoff_preseq = cutoff_preseq,
                        cutoff_chimeras = cutoff_chimeras,
                        cutoff_chrM = cutoff_chrM,
                        cutoff_cnv_mapd = cutoff_cnv_mapd,
                        cutoff_cnv_sk = cutoff_cnv_sk,
                        tight_cnv_mapd = tight_cnv_mapd,
                        tight_cnv_sk = tight_cnv_sk),
    below_floor = qc_depth_floor(df_preq$total_reads, min_value = cutoff_reads_floor),
    # CNV is not optional for this pipeline: a cell without it becomes Inconclusive rather
    # than being scored on the remaining metrics.
    required_cols = c("MAPD_CNV_Log2", "SKEW_CNV")
  )

  df_clust <- data.frame(SampleId = df_tiers$SampleId,
                         Cluster = df_tiers$TierLabel,
                         Band = df_tiers$Band,
                         stringsAsFactors = FALSE)

  #Carry the tier onto the melted CNV table so the heatmap facets by it
  melted$Cluster <- match(melted$SampleId,df_clust$SampleId) %>% df_clust$Cluster[.]

  
  p3 <- ggplot(data = melted,aes(Range,SampleId)) +
    geom_raster(aes(fill = valueFactor)) +
    facet_grid(Cluster~chr,space = "free",scales = "free") +
    theme_ohchibi(size_panel_border = 0.3) +
    scale_fill_manual(values = paleta ,name = "Estimated Ploidy") +
    theme(
      legend.position = "bottom",
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = size_labels_samples),
      axis.ticks.y = element_line(size = 0.2),
      strip.text.x = element_text(size = 9,color = "grey30"),
      strip.text.y = element_blank(),
      legend.text = element_text(size = 9,color = "grey30"),
      legend.title = element_text(size = 10),
      plot.title = element_text(size = 10),
      panel.spacing.y = unit(0.1, "lines"),
      
    ) +
    xlab(label = "Bin")  +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0))
  
  #Now start placing the graphs with the intervals
  df_preq$Cluster <- match(df_preq$SampleId,df_clust$SampleId) %>% df_clust$Cluster[.]

  # Drop unused factor levels to prevent aggregate from failing on empty levels
  df_preq$Cluster <- droplevels(df_preq$Cluster)

  # Check for empty data before aggregation
  df_preq_for_agg <- df_preq[,c("SampleId","Cluster","preseq_count")] %>% unique
  df_preq_for_agg <- df_preq_for_agg[!is.na(df_preq_for_agg$Cluster), ]
  if(nrow(df_preq_for_agg) == 0){
    stop("ERROR: No rows with valid Cluster values for aggregation.")
  }

  df_ag <- aggregate(preseq_count~Cluster, df_preq_for_agg, quantile)
  
  colnames(df_ag$preseq_count) <- colnames(df_ag$preseq_count) %>%
    gsub(pattern = "%",replacement = "") %>%
    gsub(pattern = "^",replacement = "PercPreSeq")
  
  df_preq <- match(df_preq$Cluster,df_ag$Cluster) %>% df_ag[.,-1] %>% cbind(df_preq,.)
  
  
  p_preseq <- ggplot(data = df_preq,aes(preseq_count,SampleId)) +
    geom_rect(aes(xmin =PercPreSeq25 ,xmax =PercPreSeq75,ymin = -Inf,ymax = Inf,group = Cluster),
              fill = "#A0CC2C",color = NA,alpha = 0.1)+
    geom_vline(aes(xintercept = PercPreSeq50, group = Cluster), colour = '#DD14D3',
               size = 1.3) +
    geom_vline(xintercept = cutoff_preseq,color = "red",linetype = "longdash")+
    geom_point(size = 1, alpha = 0.85)+
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
      axis.text.x = element_text(size = 7,color = "grey30",angle = 90,vjust = 0.5,hjust = 1),
      panel.background = element_blank(),
      #panel.border = element_blank(),
      axis.line = element_line(color = 'black',size = 0.3),
      strip.text.y = element_blank(),
      panel.spacing.y = unit(0.1, "lines"),
      axis.title.x = element_text(size = 9)
    ) +
    xlab(label = "Low-pass\nPreSeq value") +
    scale_x_continuous(expand = c(0,0),breaks = c(0,1000000000,2000000000,3000000000,3500000000,4000000000,5E9),
                       limits = c(0,5000000000))
  
  
  # Remove rows with NA values before aggregation
  df_ag <- df_preq[,c("SampleId","Cluster","MAPD_CNV_Log2")] %>% unique
  df_ag <- df_ag[!is.na(df_ag$MAPD_CNV_Log2), ]
  df_ag$Cluster <- droplevels(df_ag$Cluster)

  # If all values are NA, create dummy dataframe with zeros
  if(nrow(df_ag) == 0){
    warning("WARNING: No valid MAPD_CNV_Log2 values - using zeros for quantiles")
    cluster_levels <- levels(droplevels(df_preq$Cluster))
    df_ag <- data.frame(Cluster = factor(cluster_levels, levels = cluster_levels))
    df_ag$MAPD_CNV_Log2 <- matrix(rep(0, 5 * length(cluster_levels)), nrow = length(cluster_levels),
                                   dimnames = list(NULL, c("0%", "25%", "50%", "75%", "100%")))
  } else {
    df_ag <- aggregate(MAPD_CNV_Log2~Cluster,df_ag,quantile)
  }
  
  colnames(df_ag$MAPD_CNV_Log2) <- colnames(df_ag$MAPD_CNV_Log2) %>%
    gsub(pattern = "%",replacement = "") %>%
    gsub(pattern = "^",replacement = "PercCNVMAPD")
  
  df_preq <- match(df_preq$Cluster,df_ag$Cluster) %>% df_ag[.,-1] %>% cbind(df_preq,.)
  
  
  #CNV-related metrics 
  p_cnv_mapd <- ggplot(data = df_preq,aes(MAPD_CNV_Log2,SampleId)) +
    geom_rect(aes(xmin =PercCNVMAPD25 ,xmax =PercCNVMAPD75,ymin = -Inf,ymax = Inf,group = Cluster),
              fill = "#A0CC2C",color = NA,alpha = 0.1)+
    geom_vline(aes(xintercept = PercCNVMAPD50, group = Cluster), colour = '#DD14D3',
               size = 1.3) +
    geom_vline(xintercept = cutoff_cnv_mapd,color = "red",linetype = "longdash")+
    geom_point(size = 1, alpha = 0.85)+
    #geom_bar(stat = "identity",,width = 1,fill = "grey30",color = NA) +
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
      axis.text.x = element_text(size = 7,color = "grey30",angle = 90,vjust = 0.5,hjust = 1),
      panel.background = element_blank(),
      #panel.border = element_blank(),
      axis.line = element_line(color = 'black',size = 0.3),
      strip.text.y = element_blank(),
      panel.spacing.y = unit(0.1, "lines"),
      axis.title.x = element_text(size = 9)
    ) +
    xlab(label = "Low-pass\nCNV MAPD") +
    scale_x_continuous(expand = c(0,0),limits = c(0,1),oob = squish,breaks = seq(0,1,0.1))
  
  
  # Remove rows with NA values before aggregation
  df_ag <- df_preq[,c("SampleId","Cluster","SKEW_CNV")] %>% unique
  df_ag <- df_ag[!is.na(df_ag$SKEW_CNV), ]
  df_ag$Cluster <- droplevels(df_ag$Cluster)

  # If all values are NA, create dummy dataframe with zeros
  if(nrow(df_ag) == 0){
    warning("WARNING: No valid SKEW_CNV values - using zeros for quantiles")
    cluster_levels <- levels(droplevels(df_preq$Cluster))
    df_ag <- data.frame(Cluster = factor(cluster_levels, levels = cluster_levels))
    df_ag$SKEW_CNV <- matrix(rep(0, 5 * length(cluster_levels)), nrow = length(cluster_levels),
                              dimnames = list(NULL, c("0%", "25%", "50%", "75%", "100%")))
  } else {
    df_ag <- aggregate(SKEW_CNV~Cluster,df_ag,quantile)
  }
  
  colnames(df_ag$SKEW_CNV) <- colnames(df_ag$SKEW_CNV) %>%
    gsub(pattern = "%",replacement = "") %>%
    gsub(pattern = "^",replacement = "PercCNVSK")
  
  df_preq <- match(df_preq$Cluster,df_ag$Cluster) %>% df_ag[.,-1] %>% cbind(df_preq,.)
  
  
  #CNV-related metrics 
  p_cnv_skw <- ggplot(data = df_preq,aes(SKEW_CNV,SampleId)) +
    geom_rect(aes(xmin =PercCNVSK25 ,xmax =PercCNVSK75,ymin = -Inf,ymax = Inf,group = Cluster),
              fill = "#A0CC2C",color = NA,alpha = 0.1)+
    geom_vline(aes(xintercept = PercCNVSK50, group = Cluster), colour = '#DD14D3',
               size = 1.3) +
    geom_vline(xintercept = cutoff_cnv_sk,color = "red",linetype = "longdash")+
    geom_point(size = 1, alpha = 0.85)+
    #geom_bar(stat = "identity",,width = 1,fill = "grey30",color = NA) +
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
      axis.text.x = element_text(size = 7,color = "grey30",angle = 90,vjust = 0.5,hjust = 1),
      panel.background = element_blank(),
      #panel.border = element_blank(),
      axis.line = element_line(color = 'black',size = 0.3),
      strip.text.y = element_blank(),
      panel.spacing.y = unit(0.1, "lines"),
      axis.title.x = element_text(size = 9)
    ) +
    xlab(label = "Low-pass\nCNV Skewness") +
    scale_x_continuous(expand = c(0,0),limits = c(0,1),oob = squish,breaks = seq(0,1,0.1))
  
  
  #Chimeras and chrM
  
  df_ag <- df_preq[,c("SampleId","Cluster","PCT_CHIMERAS")] %>% unique %>%
    aggregate(PCT_CHIMERAS~Cluster,.,quantile)
  
  colnames(df_ag$PCT_CHIMERAS) <- colnames(df_ag$PCT_CHIMERAS) %>%
    gsub(pattern = "%",replacement = "") %>%
    gsub(pattern = "^",replacement = "PercChim")
  
  df_preq <- match(df_preq$Cluster,df_ag$Cluster) %>% df_ag[.,-1] %>% cbind(df_preq,.)
  
  
  
  p_chim <- ggplot(data = df_preq,aes(PCT_CHIMERAS,SampleId)) +
    geom_rect(aes(xmin =PercChim25 ,xmax =PercChim75,ymin = -Inf,ymax = Inf,group = Cluster),
              fill = "#A0CC2C",color = NA,alpha = 0.1)+
    geom_vline(aes(xintercept = PercChim50, group = Cluster), colour = '#DD14D3',
               size = 1.3) +
    geom_vline(xintercept = cutoff_chimeras,color = "red",linetype = "longdash")+
    geom_point(size = 1, alpha = 0.85)+
    #geom_bar(stat = "identity",,width = 1,fill = "grey30",color = NA) +
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
      axis.text.x = element_text(size = 7,color = "grey30",angle = 90,vjust = 0.5,hjust = 1),
      panel.background = element_blank(),
      #panel.border = element_blank(),
      axis.line = element_line(color = 'black',size = 0.3),
      strip.text.y = element_text(size = 9,angle = 0),
      panel.spacing.y = unit(0.1, "lines"),
      axis.title.x = element_text(size = 9)
    ) +
    xlab(label = "Low-pass\nPCT_CHIMERAS") +
    scale_x_continuous(expand = c(0,0),limits = c(0,1),oob = squish,breaks = seq(0,1,0.1))
  
  df_ag <- df_preq[,c("SampleId","Cluster","chrM")] %>% unique %>%
    aggregate(chrM~Cluster,.,quantile)
  
  colnames(df_ag$chrM) <- colnames(df_ag$chrM) %>%
    gsub(pattern = "%",replacement = "") %>%
    gsub(pattern = "^",replacement = "PercchrM")
  
  df_preq <- match(df_preq$Cluster,df_ag$Cluster) %>% df_ag[.,-1] %>% cbind(df_preq,.)
  
  
  #Check chrM
  p_chrm <- ggplot(data = df_preq,aes(chrM,SampleId)) +
    geom_rect(aes(xmin =PercchrM25 ,xmax =PercchrM75,ymin = -Inf,ymax = Inf,group = Cluster),
              fill = "#A0CC2C",color = NA,alpha = 0.1)+
    geom_vline(aes(xintercept = PercchrM50, group = Cluster), colour = '#DD14D3',
               size = 1.3) +
    geom_vline(xintercept = cutoff_chrM,color = "red",linetype = "longdash")+
    geom_point(size = 1, alpha = 0.85)+
    #geom_bar(stat = "identity",,width = 1,fill = "grey30",color = NA) +
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
      axis.text.x = element_text(size = 7,color = "grey30",angle = 90,vjust = 0.5,hjust = 1),
      panel.background = element_blank(),
      #panel.border = element_blank(),
      axis.line = element_line(color = 'black',size = 0.3),
      strip.text.y = element_blank(),
      panel.spacing.y = unit(0.1, "lines"),
      axis.title.x = element_text(size = 9)
    ) +
    xlab(label = "Low-pass\nchrM") +
    scale_x_continuous(expand = c(0,0),limits = c(0,1),oob = squish,breaks = seq(0,1,0.1))
  
  
  ### Check if metadata is passed and append it as a bar graph
  if(is.null(df_meta)){
    
  composition <- egg::ggarrange(p2,p_preseq,p_cnv_mapd,p_cnv_skw,p_chrm,p_chim,
                                    nrow =1 ,
                                    widths = c(relative_size_cnv_plot,relative_size_preseq_plot,relative_size_mapd_plot,
                                               relative_size_skew_plot,relative_size_chrm_plot,relative_size_chimeras_plot))
    dev.off()
    
    
  }else{
    
  #Create bar graph plotting group
  df_meta <- df_meta %>% 
      dplyr::filter(.data =.,SampleId %in% usable_ids) %>% droplevels
    
  df_preq$group <- match(df_preq$SampleId,df_meta$SampleId) %>% df_meta$group[.]
  df_preq$group <- factor(df_preq$group)
  df_preq$group <- droplevels(df_preq$group)
  # Dynamic palette: Tableau_20 (via paletteer) is fixed-length and fails with ggplot2 3.5+
  # when group levels exceed palette size (vec_slice OOB). Cycle brand colors.
  group_palette_base <- bskb_col[bskb_col != "#FFFFFF"]
  n_group <- nlevels(df_preq$group)
  group_fill_colors <- stats::setNames(rep_len(group_palette_base, n_group), levels(df_preq$group))

  df_preq$Bar <- "Bar"
  
  p_bar <- ggplot(data = df_preq,aes(Bar,SampleId))+
    geom_raster(aes(fill = group))+
    facet_grid(Cluster~.,space = "free",scales = "free") +
    theme_ohchibi(size_panel_border = 0.3) +
      theme(
        legend.position = "right",
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(size = 9,color = "grey30"),
        strip.text.y = element_text(size = 9,angle = 0),
        legend.text = element_text(size = 9,color = "grey30"),
        legend.title = element_text(size = 9),
        plot.title = element_text(size = 9),
        panel.spacing.y = unit(0.1, "lines"),
        
        
      )  +
      scale_x_discrete(expand = c(0,0)) +
      scale_fill_manual(values = group_fill_colors, name = "Group", drop = FALSE)
    
    
  #### Create composition
  composition <- egg::ggarrange(p3,p_preseq,p_cnv_mapd,p_cnv_skw,p_chrm,p_chim + theme(strip.text.y = element_blank()) ,p_bar,
                                nrow =1 ,
                                widths = c(relative_size_cnv_plot,relative_size_preseq_plot,relative_size_mapd_plot,
                                           relative_size_skew_plot,relative_size_chrm_plot,relative_size_chimeras_plot,relative_size_bar_plot_group))
  dev.off()
  
    
  
  }
  
  
  #Return the per-cell verdict and the composition.
  #`CompositeScore` is retained as the numeric 1-5 tier so downstream consumers and the
  #published *_ConsensusScores.txt keep a stable column name; QC_Band is the Pass/Borderline/
  #Fail collapse and BlockingMetric says which condition held each cell back.
  merged <- df_tiers %>%
    dplyr::transmute(SampleId = SampleId,
                     CompositeScore = Tier,
                     QC_Label = as.character(TierLabel),
                     QC_Band = as.character(Band),
                     BlockingMetric = BlockingMetric) %>%
    merge(df_preq %>% dplyr::select(dplyr::any_of(
            c("SampleId","preseq_count","PCT_CHIMERAS","chrM",
              "MAPD_CNV_Log2","SKEW_CNV","total_reads"))),
          by = "SampleId", all.x = TRUE) %>%
    dplyr::arrange(.data =., dplyr::desc(CompositeScore), SampleId)

  #Cells per tier, and per Pass/Borderline/Fail band
  df_sum_cat <- merged$QC_Label %>% table %>% data.frame %>%
    dplyr::rename(.data =.,Category = ".",NumberCells = Freq) %>%
    dplyr::mutate(.data =.,ProportionCells = ((NumberCells/sum(NumberCells))*100) %>% round(2))

  df_sum_band <- merged$QC_Band %>% table %>% data.frame %>%
    dplyr::rename(.data =.,Band = ".",NumberCells = Freq) %>%
    dplyr::mutate(.data =.,ProportionCells = ((NumberCells/sum(NumberCells))*100) %>% round(2))
  
  #Check if metadata was provided
  if(is.null(df_meta)){
    
    return(list(df_verdict = merged,df_sum_verdict = df_sum_cat,
                df_sum_band = df_sum_band,composition = composition))
    
    
  }else{
    
    merged$Group <- match(merged$SampleId,df_meta$SampleId) %>% df_meta$group[.]
    merged <- merged %>%
      dplyr::relocate(.data =.,c("SampleId","CompositeScore","Group"))

    df_sum_cat_group <- merged %>%
      dplyr::mutate(.data =.,Dummy =1) %>%
      reshape2::dcast(formula = Group~CompositeScore,fun.aggregate = sum,fill = 0,value.var = "Dummy")

    return(list(df_verdict = merged,
                df_sum_verdict = df_sum_cat,
                df_sum_band = df_sum_band,
                df_sum_verdict_group = df_sum_cat_group,
                composition = composition))
    
    
  }
    
  
  
}


option_list <- list(
  make_option("--seg_copy_file", type = "character", help = "SegCopy file (TSV)"),
  make_option("--metrics_file", type = "character", help = "Metrics file (TSV)"),
  make_option("--selected_metrics", type = "character", default = NULL, help = "Selected metrics file (TSV, optional)"),
  make_option("--metadata_file", type = "character", default = NULL, help = "Metadata file (TSV, optional)"),
  make_option("--cnv_summary_file", type = "character", default = NULL, help = "CNV summary file (TSV, optional, for merging)"),
  make_option("--plot_qc_config", type = "character", default = "plot_qc_config.json", help = "DNA QC config file (JSON, optional)")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$seg_copy_file) || is.null(opt$metrics_file)) {
  stop("Please provide --seg_copy_file and --metrics_file")
}

# If a CNV summary file is provided, merge it with the metrics file to create *_withcnv.txt
metrics_file_for_qc <- opt$metrics_file
if (!is.null(opt$cnv_summary_file)) {
  merged_metrics_file <- "nf-preseq-pipeline_all_metrics_mqc_withcnv.txt"
  if (!file.exists(merged_metrics_file)) {
    df_cnv <- read.table(file = opt$cnv_summary_file, header = TRUE) %>%
      dplyr::select(.data =., c("SampleId", "MAPD_CNV_Log2", "SKEW_CNV"))
    df_write <- read.table(file = opt$metrics_file, header = TRUE, sep = "\t") %>%
      dplyr::rename(.data =., SampleId = sample_name) %>%
      merge(df_cnv)
    colnames(df_write)[1] <- "sample_name"
    write.table(x = df_write, file = merged_metrics_file, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  }
  metrics_file_for_qc <- merged_metrics_file
}

# If both CNV summary file and selected metrics file are provided, merge CNV data with selected metrics
if (!is.null(opt$cnv_summary_file) && !is.null(opt$selected_metrics)) {
  merged_selected_metrics_file <- "nf-preseq-pipeline_selected_metrics_mqc_withcnv.txt"
  if (!file.exists(merged_selected_metrics_file)) {
    df_cnv <- read.table(file = opt$cnv_summary_file, header = TRUE) %>%
      dplyr::select(.data =., c("SampleId", "MAPD_CNV_Log2", "SKEW_CNV"))
    df_selected <- read.table(file = opt$selected_metrics, header = TRUE, sep = "\t") %>%
      dplyr::rename(.data =., SampleId = sample_name) %>%
      merge(df_cnv)
    colnames(df_selected)[1] <- "sample_name"
    write.table(x = df_selected, file = merged_selected_metrics_file, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  }
}

res <- plot_qc_dna(
  seg_copy_file = opt$seg_copy_file,
  metrics_file = metrics_file_for_qc,
  metadata_file = opt$metadata_file,
  # Always use config file values
  cutoff_preseq = fromJSON(opt$plot_qc_config)$cutoff_preseq,
  cutoff_chimeras = fromJSON(opt$plot_qc_config)$cutoff_chimeras,
  cutoff_chrM = fromJSON(opt$plot_qc_config)$cutoff_chrM,
  cutoff_cnv_mapd = fromJSON(opt$plot_qc_config)$cutoff_cnv_mapd,
  cutoff_cnv_sk = fromJSON(opt$plot_qc_config)$cutoff_cnv_sk,
  # tier-5 tightening and the read floor are optional in the config; fall back to the
  # function defaults so pre-existing plot_qc_config.json files keep working unchanged
  tight_cnv_mapd = fromJSON(opt$plot_qc_config)$tight_cnv_mapd %||% 0.20,
  tight_cnv_sk = fromJSON(opt$plot_qc_config)$tight_cnv_sk %||% 0.20,
  cutoff_reads_floor = fromJSON(opt$plot_qc_config)$cutoff_reads_floor %||% 1e6,
  size_labels_samples = fromJSON(opt$plot_qc_config)$size_labels_samples,
  relative_size_cnv_plot = fromJSON(opt$plot_qc_config)$relative_size_cnv_plot,
  relative_size_preseq_plot = fromJSON(opt$plot_qc_config)$relative_size_preseq_plot,
  relative_size_mapd_plot = fromJSON(opt$plot_qc_config)$relative_size_mapd_plot,
  relative_size_skew_plot = fromJSON(opt$plot_qc_config)$relative_size_skew_plot,
  relative_size_chrm_plot = fromJSON(opt$plot_qc_config)$relative_size_chrm_plot,
  relative_size_chimeras_plot = fromJSON(opt$plot_qc_config)$relative_size_chimeras_plot,
  relative_size_bar_plot_group = fromJSON(opt$plot_qc_config)$relative_size_bar_plot_group
)

write.table(res$df_verdict, "DNA-QC_ConsensusScores.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(res$df_sum_verdict, "DNA-QC_ConsensusScores_SummaryTable_mqc.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(res$df_sum_band, "DNA-QC_QCBand_SummaryTable_mqc.txt", sep = "\t", quote = FALSE, row.names = FALSE)

if (!is.null(res$df_sum_verdict_group)) {
  write.table(res$df_sum_verdict_group, "ConsensusScores_SummaryTableByGroup.txt", sep = "\t", quote = FALSE, row.names = FALSE)
}

oh.save.pdf(p = res$composition, outname = "QC_composition.pdf", outdir = "./", width = 16, height = 8)
ggsave(filename = "QC_composition_mqc.jpg", plot = res$composition, path = "./", width = 16, height = 8, units = "in", dpi = 300)

