library(dplyr)
library(magrittr)
library(data.table)

all_samples <- read.table("all_samples.txt") %$% V1

mheader <- c("",all_samples)

######################
###### NR files ######
######################

moutfile_nr <- paste0("Tab_NR.tsv")
write.table(x = mheader %>% t,file = moutfile_nr,append = FALSE,quote = FALSE,sep = "\t",row.names = FALSE,col.names = FALSE)

list_files <- read.table("list_nr.txt") %$% V1
for(mfile in list_files){

    cat("Working on ",mfile,"\n")

    df_temp <- data.table::fread(mfile,sep="\t",quote = "",data.table=FALSE)
    rownames(df_temp) <- df_temp$V1
    df_temp <- df_temp[,-1]

    Tab_ghost <- matrix(data = NA,nrow = nrow(df_temp),ncol = length(all_samples))
    rownames(Tab_ghost) <- rownames(df_temp)
    colnames(Tab_ghost) <- all_samples

    for(i in 1:ncol(Tab_ghost)){

        nom_col <- colnames(Tab_ghost)[i]
        mans <- which(colnames(df_temp) %in% nom_col)
        if(length(mans) > 0){

            Tab_ghost[,i] <- df_temp[,mans]

        }else{

            Tab_ghost[,i] <- 0

        }

    }

    write.table(x = Tab_ghost,file = moutfile_nr,append = TRUE,quote = FALSE,sep = "\t",row.names = TRUE,col.names = FALSE)
    rm(Tab_ghost)
    rm(df_temp)

}


######################
###### NR files ######
######################

moutfile_nv <- paste0("Tab_NV.tsv")
write.table(x = mheader %>% t,file = moutfile_nv,append = FALSE,quote = FALSE,sep = "\t",row.names = FALSE,col.names = FALSE)

list_files <- read.table("list_nv.txt") %$% V1
for(mfile in list_files){

    cat("Working on ",mfile,"\n")

    df_temp <- data.table::fread(mfile,sep="\t",quote = "",data.table=FALSE)
    rownames(df_temp) <- df_temp$V1
    df_temp <- df_temp[,-1]

    Tab_ghost <- matrix(data = NA,nrow = nrow(df_temp),ncol = length(all_samples))
    rownames(Tab_ghost) <- rownames(df_temp)
    colnames(Tab_ghost) <- all_samples

    for(i in 1:ncol(Tab_ghost)){

        nom_col <- colnames(Tab_ghost)[i]
        mans <- which(colnames(df_temp) %in% nom_col)
        if(length(mans) > 0){

            Tab_ghost[,i] <- df_temp[,mans]

        }else{

            Tab_ghost[,i] <- 0

        }

    }

    write.table(x = Tab_ghost,file = moutfile_nv,append = TRUE,quote = FALSE,sep = "\t",row.names = TRUE,col.names = FALSE)
    rm(Tab_ghost)
    rm(df_temp)
    
}