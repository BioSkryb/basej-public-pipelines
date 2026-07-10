library(dplyr)
library(magrittr)
library(optparse)

# Function definitions
countSegs <- function(x.vec){
  tmp.segs <- 2
  num.segs <- 0
  for(i in 1:(length(x.vec)-1)){
    if(x.vec[i]!=x.vec[i+1]){
      tmp.segs <- x.vec[i+1]
      num.segs <- num.segs+1
    }
  }
  return(num.segs)
}

rebin_segs_size <- function(df = NULL,infrabin_size = 3000000,passed_function = median){
  df$NumBases <- df[,3]-df[,2]
  mchrs <- df[,1] %>% as.character %>% unique

  # Initialize PrevNumBases column with appropriate length
  df$PrevNumBases <- NA

  # Process each chromosome separately and assign PrevNumBases in place
  for(mchr in mchrs){
    chr_indices <- which(df[,1] == mchr)
    if(length(chr_indices) > 0){
      # First bin of chromosome: use 1
      df$PrevNumBases[chr_indices[1]] <- 1
      # Subsequent bins: calculate distance from previous bin's end
      if(length(chr_indices) > 1){
        for(i in 2:length(chr_indices)){
          curr_idx <- chr_indices[i]
          prev_idx <- chr_indices[i-1]
          df$PrevNumBases[curr_idx] <- df[curr_idx, 2] - df[prev_idx, 3]
        }
      }
    }
  }
  df_new <- NULL
  for(mchr in mchrs){
    df_temp <- df %>% subset(CHR == mchr) %>% droplevels
    block_counter <- 1
    sum_block <- 0
    df_indices <- NULL
    for(i in 1:nrow(df_temp)){
      sum_block <- sum_block + df_temp$NumBases[i]
      if(sum_block < infrabin_size){
        df_indices <- data.frame(Index = i,Block = paste0("Block",block_counter)) %>%
          rbind(df_indices,.)
      }else{
        block_counter <- block_counter + 1
        df_indices <- data.frame(Index = i,Block = paste0("Block",block_counter)) %>%
          rbind(df_indices,.)
        sum_block <- df_temp$NumBases[i]
      }
    }
    mblocks <- df_temp$InfraBlock <- df_indices$Block
    vars_measure <- which(!(colnames(df_temp) %in% c("CHR","START","END","NumBases","PrevNumBases","InfraBlock")))
    for(mblock in mblocks){
      df_b <- df_temp %>% subset(InfraBlock == mblock)
      df_c <- df_b[, vars_measure, drop = FALSE] %>% apply(MARGIN =2,FUN=passed_function) %>% as.data.frame %>% t
      df_c <-data.frame(CHR = mchr,START = df_b$START %>% min,END = df_b$END %>% max) %>% cbind(df_c)
      rownames(df_c) <- NULL
      df_new <- rbind(df_new,df_c)
    }
  }
  return(df_new %>% unique)
}

# Main function
ginkgo_segment_summary <- function(rds_file, out_file) {
  Dat <- readRDS(rds_file)
  normal <- Dat$SegNorm[, 4:ncol(Dat$SegNorm), drop = FALSE]
  final <- Dat$SegCopy[, 4:ncol(Dat$SegCopy), drop = FALSE]
  size_infrabin <- 5000000
  SegNorm_rebin <- rebin_segs_size(df = Dat$SegNorm,infrabin_size = size_infrabin,passed_function = mean)   
  SegCopy_rebin <- rebin_segs_size(df = Dat$SegCopy,infrabin_size = size_infrabin,passed_function = median)   
  normal_rebin <- SegNorm_rebin[, 4:ncol(SegNorm_rebin), drop = FALSE]
  final_rebin <- SegCopy_rebin[, 4:ncol(SegCopy_rebin), drop = FALSE]
  l <- nrow(normal)
  res <- NULL
  for(k in 1:ncol(normal)){
    num_segments <- countSegs(final[,k])
    top <- 5
    mlabel <- colnames(normal)[k]
    CN <- Dat$CN[[mlabel]]
    clouds=data.frame(x=1:l, y=normal[,k]*CN)
    clouds$CN <- final[,k]
    mad_trad <- clouds$y %>% mad
    clouds$Dummy <- 1
    df_w <- aggregate(Dummy~CN,clouds,sum) %>%
      dplyr::mutate(.data =.,Prop = Dummy/sum(Dummy))
    df_w_one <- df_w %>% subset(Dummy ==1) 
    vec_one <- df_w_one$CN
    df_w_notone <- df_w %>% subset(Dummy !=1)  %>%
      dplyr::select(.data =.,c("CN","Prop"))
    df_w_one <- data.frame(CN =9999999,Prop = df_w_one$Prop %>% sum)
    df_w <- rbind(df_w_one,df_w_notone)
    df_temp_a <- clouds %>%
      dplyr::filter(.data =.,!CN %in% vec_one) %>%
      aggregate(y~CN,.,mad)
    df_temp_b <- data.frame(CN = 9999999,y = clouds %>%
                              dplyr::filter(.data =.,CN %in% vec_one) %$% y %>% mad)
    df_m <- rbind(df_temp_a,df_temp_b)
    df_m$Prop <- match(df_m$CN,df_w$CN) %>% df_w$Prop[.]
    mad_123 <- weighted.mean(df_m$y, df_m$Prop)
    sum_123 <- aggregate(Dummy~CN,clouds,sum) %>% dplyr::filter(.data=.,CN %in% c(1,2,3)) %$% Dummy %>% sum
    perc_123 <- round((sum_123/nrow(clouds))*100,3)
    clouds$DiffExpectedCN <- clouds$y - clouds$CN
    clouds$CHR <-  Dat$SegCopy[,1]
    clouds$START <-  Dat$SegCopy[,2]
    clouds$END <-  Dat$SegCopy[,3]
    vec_diff <- NULL
    vec_diff_log2 <- NULL
    vec_diff_log10 <- NULL
    for(mchr in clouds$CHR %>% unique){
      clouds_temp <- clouds %>% subset(CHR ==  mchr) %>% droplevels
      # Skip chromosomes with only 1 row (can't compute differences)
      if(nrow(clouds_temp) < 2) next
      for(i in 1:(nrow(clouds_temp)-1)){
        temp_diff <- abs((clouds_temp$y[i+1])-(clouds_temp$y[i]))
        vec_diff <- c(vec_diff,temp_diff)
        temp_diff_log2 <- abs(log2(clouds_temp$y[i+1])-log2(clouds_temp$y[i]))
        vec_diff_log2 <- c(vec_diff_log2,temp_diff_log2)
      }
    }
    clouds_rebin=data.frame(x=1:nrow(normal_rebin), y=normal_rebin[,k]*CN)
    clouds_rebin$CN <- final_rebin[,k]
    clouds_rebin$CHR <-  SegCopy_rebin[,1]
    vec_diff_rebin <- NULL
    vec_diff_log2_rebin <- NULL
    for(mchr in clouds_rebin$CHR %>% unique){
      clouds_temp <- clouds_rebin %>% subset(CHR == mchr) %>% droplevels
      # Skip chromosomes with only 1 row (can't compute differences)
      if(nrow(clouds_temp) < 2) next
      for(i in 1:(nrow(clouds_temp)-1)){
        temp_diff <- abs((clouds_temp$y[i+1])-(clouds_temp$y[i]))
        vec_diff_rebin <- c(vec_diff_rebin,temp_diff)
        temp_diff_log2 <- abs(log2(clouds_temp$y[i+1])-log2(clouds_temp$y[i]))
        vec_diff_log2_rebin <- c(vec_diff_log2_rebin,temp_diff_log2)
      }
    }
    clouds$Direction <- "Positive"
    clouds$Direction[which(clouds$DiffExpectedCN<0)] <- "Negative"
    clouds$Dummy <- 1
    df_ag <- aggregate(Dummy~CN+Direction,clouds,sum)
    mcalls <- df_ag$CN %>% table %>% data.frame %>% subset(Freq > 1) %$% . %>% as.character
    df_ratio <- NULL
    for (mcall in mcalls){
      x <- df_ag %>% subset(CN == mcall) %$% Dummy %>% sort
      df_ratio <- data.frame(W = sum(x),Ratio = x[1]/x[2]) %>%
        rbind(df_ratio,.)
    }
    df_ratio$W <- df_ratio$W/nrow(clouds)
    i_sk <- 1-(weighted.mean(df_ratio$Ratio,df_ratio$W))
    noise_naive <- mean(abs(2-clouds$y))
    noise_cn_aware <- mean(abs(clouds$DiffExpectedCN))
    noise_onlytwo <- clouds %>% subset(CN == 2) %$% DiffExpectedCN %>% abs %>% mean
    noise_naive_med <- median(abs(2-clouds$y))
    noise_cn_aware_med <- median(abs(clouds$DiffExpectedCN))
    noise_onlytwo_med <- clouds %>% subset(CN == 2) %$% DiffExpectedCN %>% abs %>% median
    res <- data.frame(SampleId = mlabel,GenomePloidy = CN,NumSegments = num_segments,SegmentMAD_Traditional = mad_trad,SegmentMAD_Aware = mad_123 ,PercDataUse123Only = perc_123,
                      MAPD_CNV = median(vec_diff),MAPD_CNV_Log2 = median(vec_diff_log2),
                      MAPD_CNV_Rebin = median(vec_diff_rebin),MAPD_CNVLog2_Rebin = median(vec_diff_log2_rebin),SKEW_CNV = i_sk,
                      Noise_Naive_Mean = noise_naive, Noise_CN_Aware_Mean = noise_cn_aware,Noise2_Mean =noise_onlytwo,
                      Noise_Naive_Median = noise_naive_med, Noise_CN_Aware_Median = noise_cn_aware_med, Noise2_Median =noise_onlytwo_med) %>%rbind(res,.)
  }
  write.table(res, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

# Command line interface
option_list <- list(
  make_option(c("--rds_file"), type = "character", default = NULL,
              help = "Path to input RDS file", metavar = "character"),
  make_option(c("--out_file"), type = "character", default = "AllSample-GinkgoSegmentSummary.txt",
              help = "Output summary file", metavar = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$rds_file)) {
  stop("Please provide --rds_file argument.")
}

ginkgo_segment_summary(opt$rds_file, opt$out_file)