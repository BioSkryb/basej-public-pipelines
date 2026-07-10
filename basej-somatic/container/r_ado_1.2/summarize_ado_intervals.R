library(ggplot2)
library(scales)


args <- commandArgs(trailingOnly=TRUE)
file_name <- args[1]
cov_cutoff <- args[2]

df_temp <- read.table(file_name,header = F,sep = "\t")
df_filtered <- df_temp[which(df_temp[,8] >= cov_cutoff),]
x <- df_filtered[,9]
a = length(x[x <0.1])
b = length(x[(x >= 0.1) &  (x < 0.2)])
c = length(x[(x >= 0.2) &  (x < 0.3)])
d = length(x[(x >= 0.3) &  (x < 0.4)])
e = length(x[(x >= 0.4) &  (x < 0.5)])
f = length(x[(x >= 0.5) &  (x < 0.6)])
g = length(x[(x >= 0.6) &  (x < 0.7)])
h = length(x[(x >= 0.7) &  (x < 0.8)])
i = length(x[(x >= 0.8) &  (x < 0.9)])
j = length(x[(x >= 0.9) &  (x <= 1.0)])
mtotal <- a+b+c+d+e+f+g+h+i+j
if (mtotal != length(x)){
	cat("Error\n")
	stop()
}

list_colnames <- c("[0-0.1)","[0.1-0.2)","[0.2-0.3)","[0.3-0.4)","[0.4-0.5)","[0.5-0.6)","[0.6-0.7)","[0.7-0.8)","[0.8-0.9)","[0.9-1]")
list_values <- c(a,b,c,d,e,f,g,h,i,j)
df_freq <- data.frame(Interval = list_colnames,Freq = list_values)
df_freq$Prop <- df_freq$Freq/mtotal
df_freq$File <- gsub(pattern ="df_ADO_",replacement = "",gsub(file_name,pattern = "\\.bqsr.*",replacement=""))

#Here change order of columns
df_freq <- df_freq[,c("File","Interval","Freq","Prop")]

write.table(df_freq,file = paste0("res_ADO_",df_freq$File[1],".tsv"),
        append = F,sep = "\t",col.names = F,row.names = F, quote =F)
  

p <- ggplot(df_freq,aes(Interval,Prop)) +
 geom_bar(stat = "identity") +
  theme_minimal() +
   ggtitle(label = file_name) +
   theme(
   axis.text.x = element_text (angle = 90,vjust = 0.5,hjust =1)
   ) +
   scale_y_continuous(breaks = seq(0,0.20,0.05),limits = c(0,0.20),oob = squish)

ggsave(plot=p,filename=paste0("histogram_ado_",df_freq$File[1],".png"),width = 6,height = 4)
