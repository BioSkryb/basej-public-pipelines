library(dplyr)
library(magrittr)


##count segs takes a stream of segment copies (SegCopy) and reports the number of segments
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


Dat <- readRDS("ginkgo_res.binsize_1000000.RDS")

normal <- Dat$SegNorm[, 4:ncol(Dat$SegNorm), drop = FALSE]
final <- Dat$SegCopy[, 4:ncol(Dat$SegCopy), drop = FALSE]

l <- nrow(normal)
res <- NULL
for(k in 1:ncol(normal)){
  
  #Calculate number of segments
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
  
  #Determine which ones have summy of 1 to pull them across
  df_w_one <- df_w %>% subset(Dummy ==1) 
  
  vec_one <- df_w_one$CN
  
  df_w_notone <- df_w %>% subset(Dummy !=1)  %>%
    dplyr::select(.data =.,c("CN","Prop"))
  
  df_w_one <- data.frame(CN =9999999,Prop = df_w_one$Prop %>% sum)

  df_w <- rbind(df_w_one,df_w_notone)
  
  #Here compute the MAD across
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
  

  res <- data.frame(SampleId = mlabel,GenomePloidy = CN,NumSegments = num_segments,SegmentMAD_Traditional = mad_trad,SegmentMAD_Aware = mad_123 ,PercDataUse123Only = perc_123 ) %>%
  rbind(res,.)
                     

}
    
write.table(res,file = "AllSample-GinkgoSegmentSummary.txt",sep = "\t",quote = FALSE,row.names = FALSE,col.names = TRUE)
