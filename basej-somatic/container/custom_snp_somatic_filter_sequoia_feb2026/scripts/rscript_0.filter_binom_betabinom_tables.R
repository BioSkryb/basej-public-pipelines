suppressMessages(library(dplyr))
suppressMessages(library(magrittr))
suppressMessages(library(data.table))

args = commandArgs(trailingOnly=TRUE)
df_file <- args[1]
binomial_threshold <- as.numeric(args[2])
beta_threshold <- as.numeric(args[3])

df <- fread(df_file,header = TRUE,sep = "\t",data.table = FALSE,quote = "")

df$Germline_qval <- df$Germline_pval %>% p.adjust(method="BH")

df$Germline_qval_log10 <- log10(df$Germline_qval)

df$Germline_filter <-  FALSE

df$Germline_filter[which(df$Germline_qval_log10 < binomial_threshold)] <- TRUE

df$Betabinomial_filter <- FALSE

df$Betabinomial_filter[which(is.na(df$Rho) | df$Rho > beta_threshold)] <- TRUE

df$Verdict <- paste0(df$Depth_filter,"_",df$Germline_filter,"_",df$Betabinomial_filter)

chosen_variants <- df %>% subset(Verdict == "TRUE_TRUE_TRUE") %$% VariantId

write.table(x =df,file = "df_verdict.txt",append = FALSE,sep = "\t",row.names=FALSE,col.names=TRUE,quote = FALSE)

write.table(x =chosen_variants,file = "chosen_variants.txt",append = FALSE,sep = "\t",row.names=FALSE,col.names=FALSE,quote = FALSE)
