library(dplyr)
library(magrittr)
library(data.table)

args = commandArgs(trailingOnly=TRUE)

threshold_as <- as.numeric(args[1])
threshold_prop <- as.numeric(args[2])
threshold_prop_bp_under <- as.numeric(args[3])
threshold_prop_bp_upper <- as.numeric(args[4])
threshold_sd_indiv <- as.numeric(args[5])
threshold_mad_indiv <- as.numeric(args[6])
threshold_sd_both <- as.numeric(args[7])
threshold_mad_both <- as.numeric(args[8])
threshold_sd_extreme <- as.numeric(args[9])
threshold_mad_extreme <- as.numeric(args[10])
# Trimmed arg chain: the group-level (Stage 2) coverage/support thresholds and the
# NumFragments thresholds were removed with their filters, so args 11-15 no longer exist.
disable_qc <- as.logical(args[11])
# disable_bppos: neutralize ONLY the positional-bias (BPPos) filter while keeping
# AS / PropClipped active. Useful for WES/capture data where read
# start positions cluster by design (probe tiling / exon edges), making the
# start-position-diversity assumption unreliable. Optional arg (backward compatible).
disable_bppos <- if (length(args) >= 12) as.logical(args[12]) else FALSE
if (is.na(disable_bppos)) disable_bppos <- FALSE
##### Sample level filtering #######

if( disable_qc == TRUE ){

    threshold_as <- 0
    threshold_prop <- 2
    threshold_prop_bp_under <- 1000
    threshold_prop_bp_upper <- 1000
    threshold_sd_indiv <- 0    
    threshold_mad_indiv <- 0
    threshold_sd_both <- 0
    threshold_mad_both <- 0
    threshold_sd_extreme <- 0
    threshold_mad_extreme <- 0

}

df <- data.table::fread(file = "df_raw_variants.tsv",sep = "\t",quote = "",header = TRUE,data.table = FALSE)

#Clean NAS and assign values
df$PROP_FRAGMENTS_BPSTART_UNDER_HQ_MQ_BQ_F[which(is.na(df$PROP_FRAGMENTS_BPSTART_UNDER_HQ_MQ_BQ_F))] <- 1
df$PROP_FRAGMENTS_BPSTART_UNDER_HQ_MQ_BQ_R[which(is.na(df$PROP_FRAGMENTS_BPSTART_UNDER_HQ_MQ_BQ_R))] <- 1
df$PROP_FRAGMENTS_BPSTART_UPPER_HQ_MQ_BQ_F[which(is.na(df$PROP_FRAGMENTS_BPSTART_UPPER_HQ_MQ_BQ_F))] <- 1
df$PROP_FRAGMENTS_BPSTART_UPPER_HQ_MQ_BQ_R[which(is.na(df$PROP_FRAGMENTS_BPSTART_UPPER_HQ_MQ_BQ_R))] <- 1
df$SD_BPSTART_FRAGMENTS_HQ_MQ_BQ_F[which(is.na(df$SD_BPSTART_FRAGMENTS_HQ_MQ_BQ_F))] <- 0
df$SD_BPSTART_FRAGMENTS_HQ_MQ_BQ_R[which(is.na(df$SD_BPSTART_FRAGMENTS_HQ_MQ_BQ_R))] <- 0
df$MAD_BPSTART_FRAGMENTS_HQ_MQ_BQ_F[which(is.na(df$MAD_BPSTART_FRAGMENTS_HQ_MQ_BQ_F))] <- 0
df$MAD_BPSTART_FRAGMENTS_HQ_MQ_BQ_R[which(is.na(df$MAD_BPSTART_FRAGMENTS_HQ_MQ_BQ_R))] <- 0


df_as <- df %>% subset(MEDIAN_AS_VARIANT_READS >= threshold_as) %>% 
    droplevels %>% dplyr::mutate(.data =.,Verdict = "PASSED_AS") %>%
    dplyr::select(.data =.,c("Verdict","SampleId","VariantId"))

df_prop <- df %>% subset(PROP_BASES_CLIPPED < threshold_prop) %>% 
    droplevels %>% dplyr::mutate(.data =.,Verdict = "PASSED_PROPCLIPPED") %>%
    dplyr::select(.data =.,c("Verdict","SampleId","VariantId"))


#Subset by positions

df_a <- df %>% 
    subset(NUM_FRAGMENTS_HQ_MQ_BQ_F  <2 & NUM_FRAGMENTS_HQ_MQ_BQ_R > 1) %>%
    subset( (PROP_FRAGMENTS_BPSTART_UNDER_HQ_MQ_BQ_R <  threshold_prop_bp_under &  PROP_FRAGMENTS_BPSTART_UPPER_HQ_MQ_BQ_R < threshold_prop_bp_upper ) | (SD_BPSTART_FRAGMENTS_HQ_MQ_BQ_R > threshold_sd_indiv & MAD_BPSTART_FRAGMENTS_HQ_MQ_BQ_R > threshold_mad_indiv ) ) %>%
droplevels

df_b <- df %>% 
    subset(NUM_FRAGMENTS_HQ_MQ_BQ_F  > 1 & NUM_FRAGMENTS_HQ_MQ_BQ_R < 2) %>%
    subset( (PROP_FRAGMENTS_BPSTART_UNDER_HQ_MQ_BQ_F <  threshold_prop_bp_under &  PROP_FRAGMENTS_BPSTART_UPPER_HQ_MQ_BQ_F < threshold_prop_bp_upper ) | (SD_BPSTART_FRAGMENTS_HQ_MQ_BQ_F > threshold_sd_indiv & MAD_BPSTART_FRAGMENTS_HQ_MQ_BQ_F > threshold_mad_indiv ) ) %>%
droplevels

temp_c <- df %>%
    subset(NUM_FRAGMENTS_HQ_MQ_BQ_F  > 1 & NUM_FRAGMENTS_HQ_MQ_BQ_R > 1) %>% droplevels

df_c_1 <- temp_c %>%
    subset( (PROP_FRAGMENTS_BPSTART_UNDER_HQ_MQ_BQ_F < threshold_prop_bp_under & PROP_FRAGMENTS_BPSTART_UPPER_HQ_MQ_BQ_F <  threshold_prop_bp_upper) & ( (SD_BPSTART_FRAGMENTS_HQ_MQ_BQ_F > threshold_sd_both & MAD_BPSTART_FRAGMENTS_HQ_MQ_BQ_F > threshold_mad_both ) | (SD_BPSTART_FRAGMENTS_HQ_MQ_BQ_R > threshold_sd_extreme & MAD_BPSTART_FRAGMENTS_HQ_MQ_BQ_R > threshold_mad_extreme  )  ) ) %>%
    droplevels

df_c_2 <- temp_c %>%
    subset( (PROP_FRAGMENTS_BPSTART_UNDER_HQ_MQ_BQ_R < threshold_prop_bp_under & PROP_FRAGMENTS_BPSTART_UPPER_HQ_MQ_BQ_R <  threshold_prop_bp_upper) & ( (SD_BPSTART_FRAGMENTS_HQ_MQ_BQ_R > threshold_sd_both & MAD_BPSTART_FRAGMENTS_HQ_MQ_BQ_R > threshold_mad_both ) | (SD_BPSTART_FRAGMENTS_HQ_MQ_BQ_F > threshold_sd_extreme & MAD_BPSTART_FRAGMENTS_HQ_MQ_BQ_F > threshold_mad_extreme  )  ) ) %>%
    droplevels

if( disable_qc == TRUE || disable_bppos == TRUE ){

    # BPPos disabled: every cell passes the positional-bias gate (pass-through),
    # so the Dummy==3 AND reduces to AS & PropClipped.
    df_passed_BPPOS <- df %>%
        dplyr::mutate(.data =.,Verdict = "PASSED_BPPOSINREADS") %>%
        dplyr::select(.data =.,c("Verdict","SampleId","VariantId")) %>%
        unique

}else{


    df_passed_BPPOS <- rbind(df_a,df_b) %>%
        rbind(df_c_1) %>%
        rbind(df_c_2) %>%
        dplyr::mutate(.data =.,Verdict = "PASSED_BPPOSINREADS") %>%
        dplyr::select(.data =.,c("Verdict","SampleId","VariantId")) %>%
        unique

}

# NumFragments filter removed: the per-cell artifact QC is now AS & PropClipped & BPPos only.
df_end <- rbind(df_as,df_prop,df_passed_BPPOS)

#Here evaluate if we have enough variants to pass the filters
# A variant is kept if >=1 cell passed all THREE artifact filters (Dummy == 3).
if(nrow(df_end)>0){
    df_end$Dummy <- 1
    df_suma <- aggregate(Dummy~VariantId+SampleId,df_end,sum)
    variants_to_use_group <- df_suma %>% subset(Dummy == 3) %$% VariantId %>% as.character %>% unique
    num_samples <- df$SampleId %>% as.character %>% unique %>% length
}else{
    variants_to_use_group <- NULL
    num_samples <- 0
}


##### Group level filtering #######
if(length(variants_to_use_group)>0){

    #Get all the alt ids
    df_raw_variants <-  data.table::fread(file = "res_end.tsv",sep = "\t",quote = "",header = TRUE,data.table = FALSE)

    df_raw_variants_group_alt <- df_raw_variants %>%
        dplyr::filter(.data =.,VariantId %in% variants_to_use_group) %>% droplevels

    pos_to_use_group <- variants_to_use_group %>% strsplit(split = "_") %>% lapply(FUN = function(x)x[2]) %>% unlist %>% unique 

    df_raw_variants_group_ref <- df_raw_variants %>% 
        subset(ALT == "REF") %>%
        dplyr::filter(.data =.,POS %in% pos_to_use_group) %>% droplevels

    df_raw_variants_group <- rbind(df_raw_variants_group_ref,df_raw_variants_group_alt) %>% unique

    # Group-level (Stage 2) filtering removed: no position-coverage gate and no
    # variant-support gate. Every variant that passed the per-cell artifact filters
    # (AS & PropClipped & BPPos) in >=1 cell is kept. df_raw_variants_group is still
    # built above because the NR (depth) matrix is constructed from it below.
    variants_to_use_passed <- variants_to_use_group
}else{
    variants_to_use_passed <- NULL
}

if(length(variants_to_use_passed)>0){

    #Subset the chosen variants 
    df_nv <- df_raw_variants %>% 
        dplyr::filter(.data =.,VariantId %in% variants_to_use_passed) %>% droplevels %>%
        dplyr::mutate(.data =.,NUM_FRAGMENTS_HQ_MQ_BQ_FR= NUM_FRAGMENTS_HQ_MQ_BQ_F + NUM_FRAGMENTS_HQ_MQ_BQ_R) %>%
        dplyr::select(.data =.,c("SampleId","POS","VariantId","NUM_FRAGMENTS_HQ_MQ_BQ_FR")) %>% droplevels
        
    pos_nv <- df_nv$POS %>% unique

    df_nr <- df_raw_variants_group %>%
        dplyr::select(.data =.,c("SampleId","CHROM","POS","NUM_FRAGMENTS_HQ_POSITION")) %>% unique %>%
        dplyr::mutate(CHROM_POS = paste0(CHROM, "_", POS)) %>%
        dplyr::filter(.data =.,POS %in% pos_nv) %>% droplevels

    ##### Matrix creation ##############

    Tab_nv <- df_nv %>% reshape2::acast(VariantId~SampleId,value.var = "NUM_FRAGMENTS_HQ_MQ_BQ_FR",fill = 0)
    # acast returns a vector when there is only one VariantId; keep as 1-row matrix
    if (!is.matrix(Tab_nv)) {
        Tab_nv <- as.matrix(t(Tab_nv))
        rownames(Tab_nv) <- as.character(unique(df_nv$VariantId))
    }

    Tab_nr <- df_nr %>% reshape2::acast(CHROM_POS~SampleId,value.var = "NUM_FRAGMENTS_HQ_POSITION",fill = 0)
    # acast returns a vector when there is only one CHROM_POS; keep as 1-row matrix
    if (!is.matrix(Tab_nr)) {
        Tab_nr <- as.matrix(t(Tab_nr))
        rownames(Tab_nr) <- as.character(unique(df_nr$CHROM_POS))
    }

    nv_names <- rownames(Tab_nv)
    nv_chrom_pos <- sapply(strsplit(nv_names, "_"), function(x) paste0(x[1], "_", x[2]))

    # Tab_nr row names are already CHROM_POS — match directly.
    idx_nr <- match(nv_chrom_pos, rownames(Tab_nr))
    Tab_nr <- Tab_nr[idx_nr, , drop = FALSE]

    rownames(Tab_nr) <- rownames(Tab_nv)

    # Align columns: expand NV with zeros for non-carrier samples that have REF-row
    # coverage in NR, then order both matrices to the same full sample set.
    all_samples <- union(colnames(Tab_nv), colnames(Tab_nr))
    missing_in_nv <- setdiff(all_samples, colnames(Tab_nv))
    if (length(missing_in_nv) > 0) {
        zero_cols <- matrix(0L, nrow = nrow(Tab_nv), ncol = length(missing_in_nv),
                            dimnames = list(rownames(Tab_nv), missing_in_nv))
        Tab_nv <- cbind(Tab_nv, zero_cols)
    }
    Tab_nv <- Tab_nv[, all_samples, drop = FALSE]
    Tab_nr <- Tab_nr[, all_samples, drop = FALSE]


    moutfile_nr <- paste0("Tab_NR.tsv")

    moutfile_nv <- paste0("Tab_NV.tsv")

    #Write matrices
    mheader <- c("",colnames(Tab_nr))
        
    write.table(x = mheader %>% t,file = moutfile_nr,append = FALSE,quote = FALSE,sep = "\t",row.names = FALSE,col.names = FALSE)
    
    mheader <- c("",colnames(Tab_nv))

    write.table(x = mheader %>% t,file = moutfile_nv,append = FALSE,quote = FALSE,sep = "\t",row.names = FALSE,col.names = FALSE)

    write.table(x = Tab_nr,file = moutfile_nr,append = TRUE,quote = FALSE,sep = "\t",row.names = TRUE,col.names = FALSE)
        
    write.table(x = Tab_nv,file = moutfile_nv,append = TRUE,quote = FALSE,sep = "\t",row.names = TRUE,col.names = FALSE)


}else{


    mheader <- c("",df$SampleId %>% as.character %>% unique )


    moutfile_nr <- paste0("Tab_NR.tsv")

    moutfile_nv <- paste0("Tab_NV.tsv")

    #Here only write the column names
    write.table(x = mheader %>% t,file = moutfile_nr,append = FALSE,quote = FALSE,sep = "\t",row.names = FALSE,col.names = FALSE)
    write.table(x = mheader %>% t,file = moutfile_nv,append = FALSE,quote = FALSE,sep = "\t",row.names = FALSE,col.names = FALSE)
}

rm(df)
gc()

# --- All-variants matrices (no QC filtering) ---
df_all <- data.table::fread(file = "res_end.tsv", sep = "\t", quote = "", header = TRUE, data.table = FALSE)

all_variant_ids <- df_all %>%
    dplyr::filter(ALT != "REF") %$% VariantId %>% as.character %>% unique

if (length(all_variant_ids) > 0) {

    df_nv_all <- df_all %>%
        dplyr::filter(VariantId %in% all_variant_ids) %>%
        dplyr::mutate(NUM_FRAGMENTS_HQ_MQ_BQ_FR = NUM_FRAGMENTS_HQ_MQ_BQ_F + NUM_FRAGMENTS_HQ_MQ_BQ_R) %>%
        dplyr::select(SampleId, POS, VariantId, NUM_FRAGMENTS_HQ_MQ_BQ_FR)

    pos_all <- df_nv_all$POS %>% unique

    df_nr_all <- df_all %>%
        dplyr::select(SampleId, CHROM, POS, NUM_FRAGMENTS_HQ_POSITION) %>% unique %>%
        dplyr::mutate(CHROM_POS = paste0(CHROM, "_", POS)) %>%
        dplyr::filter(POS %in% pos_all)

    Tab_nv_all <- df_nv_all %>%
        reshape2::acast(VariantId ~ SampleId, value.var = "NUM_FRAGMENTS_HQ_MQ_BQ_FR", fill = 0)
    if (!is.matrix(Tab_nv_all)) {
        Tab_nv_all <- as.matrix(t(Tab_nv_all))
        rownames(Tab_nv_all) <- as.character(unique(df_nv_all$VariantId))
    }

    Tab_nr_all <- df_nr_all %>%
        reshape2::acast(CHROM_POS ~ SampleId, value.var = "NUM_FRAGMENTS_HQ_POSITION", fill = 0)
    if (!is.matrix(Tab_nr_all)) {
        Tab_nr_all <- as.matrix(t(Tab_nr_all))
        rownames(Tab_nr_all) <- as.character(unique(df_nr_all$CHROM_POS))
    }

    nv_names_all <- rownames(Tab_nv_all)
    nv_chrom_pos_all <- sapply(strsplit(nv_names_all, "_"), function(x) paste0(x[1], "_", x[2]))
    # Tab_nr_all row names are already CHROM_POS — match directly.
    idx_nr_all <- match(nv_chrom_pos_all, rownames(Tab_nr_all))
    Tab_nr_all <- Tab_nr_all[idx_nr_all, , drop = FALSE]
    rownames(Tab_nr_all) <- rownames(Tab_nv_all)

    all_samples_all <- union(colnames(Tab_nv_all), colnames(Tab_nr_all))
    missing_in_nv_all <- setdiff(all_samples_all, colnames(Tab_nv_all))
    if (length(missing_in_nv_all) > 0) {
        zero_cols <- matrix(0L, nrow = nrow(Tab_nv_all), ncol = length(missing_in_nv_all),
                            dimnames = list(rownames(Tab_nv_all), missing_in_nv_all))
        Tab_nv_all <- cbind(Tab_nv_all, zero_cols)
    }
    Tab_nv_all <- Tab_nv_all[, all_samples_all, drop = FALSE]
    Tab_nr_all <- Tab_nr_all[, all_samples_all, drop = FALSE]

    write.table(x = t(c("", colnames(Tab_nv_all))), file = "Tab_NV_all.tsv", append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    write.table(x = Tab_nv_all, file = "Tab_NV_all.tsv", append = TRUE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = FALSE)

    write.table(x = t(c("", colnames(Tab_nr_all))), file = "Tab_NR_all.tsv", append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    write.table(x = Tab_nr_all, file = "Tab_NR_all.tsv", append = TRUE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = FALSE)

} else {

    all_samples_header <- df_all$SampleId %>% as.character %>% unique
    write.table(x = t(c("", all_samples_header)), file = "Tab_NV_all.tsv", append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    write.table(x = t(c("", all_samples_header)), file = "Tab_NR_all.tsv", append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

}

#Write the quality files without header
write.table(x = df_as,file = "df_passed_AS.tsv",append = FALSE,quote = FALSE,sep = "\t",row.names = FALSE,col.names = FALSE)

write.table(x = df_prop,file = "df_passed_propclipped.tsv",append = FALSE,quote = FALSE,sep = "\t",row.names = FALSE,col.names = FALSE)

write.table(x = df_passed_BPPOS,file = "df_passed_BPPOS.tsv",append = FALSE,quote = FALSE,sep = "\t",row.names = FALSE,col.names = FALSE)
