library(dplyr)
library(magrittr)
args <- commandArgs(trailingOnly=TRUE)
rds_file <- args[1]

Dat <- readRDS(rds_file)

#Create merge structure in which we process all datasets
mids <- names(Dat$CNV_structures)


rectangles1 <- NULL
rectangles2 <- NULL
clouds <- data.frame(x = NA,y = NA,SampleId = NA)
flat <- data.frame(x = NA,y = NA,SampleId = NA)
amp <- data.frame(x = NA,y = NA,SampleId = NA)
del <- data.frame(x = NA,y = NA,SampleId = NA)
anno <- data.frame(x = NA,y = NA,chrom = NA,SampleId = NA)




for(mid in mids){
  
  Dat_sub <- Dat$CNV_structures[[mid]]
  
  rectangles1 <- Dat_sub$rectangles1 %>%
    dplyr::mutate(.data = .,SampleId = mid %>% gsub(pattern = "_sorted.*",replacement = "")) %>%
    rbind(rectangles1,.)
  
  rectangles2 <- Dat_sub$rectangles2 %>%
    dplyr::mutate(.data = .,SampleId = mid %>% gsub(pattern = "_sorted.*",replacement = "")) %>%
    rbind(rectangles2,.)
  
  clouds <- Dat_sub$clouds %>%
    dplyr::mutate(.data = .,SampleId = mid %>% gsub(pattern = "_sorted.*",replacement = "")) %>%
    rbind(clouds,.)
  
  flat <- Dat_sub$flat %>%
    dplyr::mutate(.data = .,SampleId = mid %>% gsub(pattern = "_sorted.*",replacement = "")) %>%
    rbind(flat,.)
  
  amp <- Dat_sub$amp %>%
    dplyr::mutate(.data = .,SampleId = mid %>% gsub(pattern = "_sorted.*",replacement = "")) %>%
    rbind(amp,.)
  
  del <- Dat_sub$del %>%
    dplyr::mutate(.data = .,SampleId = mid %>% gsub(pattern = "_sorted.*",replacement = "")) %>%
    rbind(del,.)
  
  anno <- Dat_sub$anno %>%
    dplyr::mutate(.data = .,SampleId = mid %>% gsub(pattern = "_sorted.*",replacement = "")) %>%
    rbind(anno,.)
}


mchrs <- anno$chrom %>% na.omit %>% unique %>% paste0("chr",.)    

num_r1 <- nrow(rectangles1)/(rectangles1$SampleId %>% unique %>% length)
num_r2 <- nrow(rectangles2)/(rectangles2$SampleId %>% unique %>% length)

flag <- 0
i <- -1
j <- 0
chr_r1 <- NULL
while (flag == 0){
	 i <- i +2
 chr_r1 <- c(chr_r1,mchrs[i])
  j <- j +1
  if(j == num_r1){
	      flag <- 1
   }
}

chr_r2 <- mchrs[which(!(mchrs %in% chr_r1))]

rectangles1$chr <- chr_r1
rectangles2$chr <- chr_r2


rectangles <- rbind(rectangles1,rectangles2)

rectangles$chr <- rectangles$chr %>%
  gsub(pattern = "chr",replacement = "")

rectangles$chr <- rectangles$chr %>% 
  factor(levels = mchrs %>% gsub(pattern = "chr",replacement = ""))

rectangles <- with(rectangles,order(chr)) %>%
  rectangles[.,]

colnames(rectangles)[1:2] <- c("start","end")

clouds$chr <- NA
flat$chr <- NA
amp$chr <- NA
del$chr <- NA


for(i in 1:nrow(rectangles)){
  
  mstart <- rectangles[i,] %$% start %>% as.numeric
  mend <- rectangles[i,] %$% end %>% as.numeric
  mchr <- rectangles[i,] %$% chr %>% as.character
  
  clouds$chr[which(dplyr::between(x = clouds$x,left = mstart,right = mend))] <- mchr
  flat$chr[which(dplyr::between(x = flat$x,left = mstart,right = mend))] <- mchr
  amp$chr[which(dplyr::between(x = amp$x,left = mstart,right = mend))] <- mchr
  del$chr[which(dplyr::between(x = del$x,left = mstart,right = mend))] <- mchr
  
}

clouds$chr <- clouds$chr %>%
  factor(levels = rectangles$chr %>% levels)

flat$chr <- flat$chr %>%
  factor(levels = rectangles$chr %>% levels)

amp$chr <- amp$chr %>%
  factor(levels = rectangles$chr %>% levels)

del$chr <- del$chr %>%
  factor(levels = rectangles$chr %>% levels)
  
#Remove NA dummy levels
clouds <- clouds %>% dplyr::filter(.data =.,!is.na(SampleId ))
flat <- flat %>% dplyr::filter(.data =.,!is.na(SampleId ))
amp <- amp %>% dplyr::filter(.data =.,!is.na(SampleId ))
del <- del %>% dplyr::filter(.data =.,!is.na(SampleId ))
anno <- anno %>% dplyr::filter(.data =.,!is.na(SampleId ))


#Write tsv files
write.table(rectangles,file = "rectangles.tsv",sep = "\t",quote = FALSE,row.names = FALSE,col.names = TRUE,append = FALSE)
write.table(clouds,file = "clouds.tsv",sep = "\t",quote = FALSE,row.names = FALSE,col.names = TRUE,append = FALSE)
write.table(flat,file = "flat.tsv",sep = "\t",quote = FALSE,row.names = FALSE,col.names = TRUE,append = FALSE)
write.table(amp,file = "amp.tsv",sep = "\t",quote = FALSE,row.names = FALSE,col.names = TRUE,append = FALSE)
write.table(del,file = "del.tsv",sep = "\t",quote = FALSE,row.names = FALSE,col.names = TRUE,append = FALSE)

#Process overall extimated sample ploidy 
df_sample_cn <- NULL
for(mid in mids){
    mcn <- Dat$CN[[mid]]
    df_sample_cn <- data.frame(SampleId = mid,SamplePloidy = mcn) %>%
      rbind(df_sample_cn,.)
}

write.table(df_sample_cn,file = "sampleploidy.tsv",sep = "\t",quote = FALSE,row.names = FALSE,col.names = TRUE,append = FALSE)
