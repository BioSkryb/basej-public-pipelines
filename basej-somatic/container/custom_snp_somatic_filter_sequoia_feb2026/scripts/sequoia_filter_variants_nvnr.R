
if(!require("optparse", character.only=T,quietly = T, warn.conflicts = F)){
  install.packages("optparse",repos = "http://cran.us.r-project.org")
  library("optparse", character.only=T,quietly = T, warn.conflicts = F)
}
#----------------------------------
# Input options
#----------------------------------
option_list = list(
  make_option(c("-v", "--input_nv"), action="store", default=NULL, type='character', help="Input NV matrix (rows are variants, columns are samples)"),
  make_option(c("-r", "--input_nr"), action="store", default=NULL, type='character', help="Input NR matrix (rows are variants, columns are samples)"),
  make_option(c("-o", "--output_dir"), action="store", default="", type='character', help="Output directory for files"),
  make_option(c("-n", "--ncores"), action="store", default=1, type='numeric', help="Number of cores to use for the beta-binomial step"),
  make_option(c("--min_cov"), action="store", default=10, type='numeric', help="Lower threshold for mean coverage across variant site"),
  make_option(c("--max_cov"), action="store", default=500, type='numeric', help="Upper threshold for mean coverage across variant site"),
  make_option(c("--cnv_samples"), action="store", default=NULL, type='character', help="Samples with CNVs, exclude from germline/depth-based filtering, separate with a comma"),
  make_option(c("--gender_passed"), action="store", default=NULL, type='character', help="Either male or female. Helps improve performance by splitting")

)
opt = parse_args(OptionParser(option_list=option_list, add_help_option=T))

print(opt)

ncores=opt$ncores
min_cov=opt$min_cov
max_cov=opt$max_cov
output_dir=opt$output_dir
if(is.null(opt$cnv_samples)) {samples_with_CNVs=NULL} else {samples_with_CNVs=unlist(strsplit(x=opt$cnv_samples,split = ","))}
nv_path=opt$input_nv
nr_path=opt$input_nr
gender_passed=opt$gender_passed

#Do not filter out anything so define rho to 1
snv_rho <- 1
indel_rho <- 1 
patient_ID <- "df_binomial_betabinomial"

#----------------------------------
# Load packages (install if they are not installed yet)
#----------------------------------
options(stringsAsFactors = F)
cran_packages=c("ggplot2","ape","seqinr","stringr","data.table","tidyr","dplyr","VGAM","MASS","devtools","parallel")
bioconductor_packages=c("Rsamtools","GenomicRanges")

for(package in cran_packages){
  if(!require(package, character.only=T,quietly = T, warn.conflicts = F)){
    install.packages(as.character(package),repos = "http://cran.us.r-project.org")
    library(package, character.only=T,quietly = T, warn.conflicts = F)
  }
}
if (!require("BiocManager", quietly = T, warn.conflicts = F))
  install.packages("BiocManager")
for(package in bioconductor_packages){
  if(!require(package, character.only=T,quietly = T, warn.conflicts = F)){
    BiocManager::install(as.character(package))
    library(package, character.only=T,quietly = T, warn.conflicts = F)
  }
}
if(!require("treemut", character.only=T,quietly = T, warn.conflicts = F)){
  install_git("https://github.com/NickWilliamsSanger/treemut")
  library("treemut",character.only=T,quietly = T, warn.conflicts = F)
}

#----------------------------------
# Functions
#----------------------------------

exact.binomial=function(gender,NV,NR,cutoff=-5,qval_return=F){
  # Function to filter out germline variants based on unmatched
  # variant calls of multiple samples from same individual (aggregate coverage
  # ideally >150 or so, but will work with less). NV is matrix of reads supporting 
  # variants and NR the matrix with total depth (samples as columns, mutations rows, 
  # with rownames as chr_pos_ref_alt or equivalent). Returns a logical vector, 
  # TRUE if mutation is likely to be germline.
  
  #Check if it is auto
  XY_chromosomal = grepl("X|Y",rownames(NR))
  flag_is_xy <- FALSE
  #Dealing with XY if higher than 0
  if(length(which(XY_chromosomal)==TRUE)>0){

    flag_is_xy <- TRUE
    XY_chromosomal <- rep(TRUE,nrow(NR))

  }else{

    autosomal = rep(TRUE,nrow(NR))
  }
  
  if(gender=="female"){

    cat("Entered female binomial filtering. FLAG_XY: ",flag_is_xy,"\n")

    NV_vec = rowSums(NV)
    NR_vec = rowSums(NR)
    pval = rep(1,length(NV_vec))
    for (n in 1:length(NV_vec)){
      if(NR_vec[n]>0){
        pval[n] = binom.test(x=NV_vec[n],
                             n=NR_vec[n],
                             p=0.5,alt='less')$p.value
      }
    }
  }
  # For male, split test in autosomal and XY chromosomal part
  if(gender=="male"){

    pval=rep(1,nrow(NV))

    if(flag_is_xy== FALSE){

      cat("Entered male autosomal binomial filtering. FLAG_XY: ",flag_is_xy,"\n")

      pval=rep(1,nrow(NV))
      NV_vec = rowSums(NV)[autosomal]
      NR_vec = rowSums(NR)[autosomal]
      pval_auto = rep(1,sum(autosomal))

      for (n in 1:sum(autosomal)){
        if(NR_vec[n]>0){
          pval_auto[n] = binom.test(x=NV_vec[n],
                                    n=NR_vec[n],
                                    p=0.5,alt='less')$p.value
        }
      }

      pval[autosomal]=pval_auto

    }else{
        
        cat("Entered male XY chromosomal binomial filtering. FLAG_XY: ",flag_is_xy,"\n")

        pval=rep(1,nrow(NV))
        
        pval_XY = rep(1,sum(XY_chromosomal))

        NV_vec = rowSums(NV)[XY_chromosomal]
        NR_vec = rowSums(NR)[XY_chromosomal]

        for (n in 1:sum(XY_chromosomal)){
          if(NR_vec[n]>0){
            pval_XY[n] = binom.test(x=NV_vec[n],
                                    n=NR_vec[n],
                                    p=0.95,alt='less')$p.value
          }
        }

        pval[XY_chromosomal]=pval_XY
    }

  }
  qval = p.adjust(pval,method="none")
  if(qval_return){
    return(qval)
  }else{
    germline = log10(qval)>cutoff
    return(germline)
  }
}

estimateRho_gridml = function(NV_vec,NR_vec) {
  # Function to estimate maximum likelihood value of rho for beta-binomial
  rhovec = 10^seq(-6,-0.05,by=0.05) # rho will be bounded within 1e-6 and 0.89
  mu=sum(NV_vec)/sum(NR_vec)
  ll = sapply(rhovec, function(rhoj) sum(dbetabinom(x=NV_vec, size=NR_vec, rho=rhoj, prob=mu, log=T)))
  return(rhovec[ll==max(ll)][1])
}

beta.binom.filter = function(NR,NV){
  # Function to apply beta-binomial filter for artefacts. Works best on sets of
  # clonal samples (ideally >10 or so). As before, takes NV and NR as input. 
  # Optionally calculates pvalue of likelihood beta-binomial with estimated rho
  # fits better than binomial. This was supposed to protect against low-depth variants,
  # but use with caution. Returns logical vector with good variants = TRUE
  
  rho_est = pval = rep(NA,nrow(NR))
  for (k in 1:nrow(NR)){
    rho_est[k]=estimateRho_gridml(NV_vec = as.numeric(NV[k,]),
                                  NR_vec=as.numeric(NR[k,]))
  }
  return(rho_est)
}

dbinomtrunc = function(x, size, prob, minx=4) {
  dbinom(x, size, prob) / pbinom(minx-0.1, size, prob, lower.tail=F)
}

estep = function(x,size,p.vector,prop.vector,ncomp, mode){
  ## p.vector = vector of probabilities for the individual components
  ## prop.vector = vector of proportions for the individual components
  ## ncomp = number of components
  p.mat_estep = matrix(0,ncol=ncomp,nrow=length(x))
  for (i in 1:ncomp){
    if(mode=="Truncated") p.mat_estep[,i]=prop.vector[i]*dbinomtrunc(x,size,prob=p.vector[i])
    if(mode=="Full") p.mat_estep[,i]=prop.vector[i]*dbinom(x,size,prob=p.vector[i])
  }
  norm = rowSums(p.mat_estep) ## normalise the probabilities
  p.mat_estep = p.mat_estep/norm
  LL = sum(log(norm)) ## log-likelihood
  
  ## classification of observations to specific components (too crude?)
  which_clust = rep(1,length(x))
  if(ncomp>1){
    which_clust = apply(p.mat_estep, 1, which.max)
  }
  
  list("posterior"=p.mat_estep,
       "LL"=LL,
       "Which_cluster"=which_clust)
}

mstep = function(x,size,e.step){
  # estimate proportions
  prop.vector_temp = colMeans(e.step$posterior)
  # estimate probabilities
  p.vector_temp = colSums(x/size*e.step$posterior) / colSums(e.step$posterior)
  
  list("prop"=prop.vector_temp,
       "p"=p.vector_temp)   
}

em.algo = function(x,size,prop.vector_inits,p.vector_inits,maxit=5000,tol=1e-6,nclust,binom_mode){
  ## prop.vector_inits =  initial values for the mixture proportions
  ## p.vector_inits =  initial values for the probabilities 
  
  # Initiate EM
  flag = 0
  e.step = estep(x,size,p.vector = p.vector_inits,prop.vector = prop.vector_inits,ncomp=nclust,mode=binom_mode)
  m.step = mstep(x,size,e.step)
  prop_cur = m.step[["prop"]]
  p_cur = m.step[["p"]]
  cur.LL = e.step[["LL"]]
  LL.vector = e.step[["LL"]]
  
  # Iterate between expectation and maximisation steps
  for (i in 2:maxit){
    e.step = estep(x,size,p.vector = p_cur,prop.vector = prop_cur,ncomp=nclust,mode=binom_mode)
    m.step = mstep(x,size,e.step)
    prop_new = m.step[["prop"]]
    p_new = m.step[["p"]]
    
    LL.vector = c(LL.vector,e.step[["LL"]])
    LL.diff = abs((cur.LL - e.step[["LL"]]))
    which_clust = e.step[["Which_cluster"]]
    # Stop iteration if the difference between the current and new log-likelihood is less than a tolerance level
    if(LL.diff < tol){ flag = 1; break}
    
    # Otherwise continue iteration
    prop_cur = prop_new; p_cur = p_new; cur.LL = e.step[["LL"]]
    
  }
  if(!flag) warning("Didn’t converge\n")
  
  BIC = log(length(x))*nclust*2-2*cur.LL
  AIC = 4*nclust-2*cur.LL
  list("LL"=LL.vector,
       "prop"=prop_cur,
       "p"=p_cur,
       "BIC"=BIC,
       "AIC"=AIC,
       "n"=nclust,
       "Which_cluster"=which_clust)
}

binom_mix = function(x,size,nrange=1:3,criterion="BIC",maxit=5000,tol=1e-6, mode="Full"){
  ## Perform the EM algorithm for different numbers of components
  ## Select best fit using the Bayesian Information Criterion (BIC) 
  ## or the Akaike information criterion (AIC)
  i=1
  results = list()
  BIC_vec = c()
  AIC_vec = c()
  
  for (n in nrange){
    ## Initialise EM algorithm with values from kmeans clustering
    init = kmeans(x/size,n)
    prop_init = init$size/length(x)
    p_init = init$centers
    
    results[[i]] = em.algo(x,size,prop.vector_inits = prop_init,p.vector_inits=p_init,nclust=n,maxit,tol,binom_mode=mode)
    BIC_vec = c(BIC_vec,results[[i]]$BIC)
    AIC_vec = c(AIC_vec,results[[i]]$AIC)
    i=i+1
  }
  if (criterion=="BIC"){
    results[[which.min(BIC_vec)]]$BIC_vec=BIC_vec
    return(results[[which.min(BIC_vec)]])
  }
  if (criterion=="AIC"){
    return(results[[which.min(AIC_vec)]])
  }
}

#----------------------------------
# Read in data
#----------------------------------
print("Reading in data...")
  
if(!is.null(nr_path)&!is.null(nv_path)){
  NR = fread(nr_path,data.table=F)
  rownames(NR)=NR[,1]
  NR=NR[,-1]
  samples_exclude_int <- names(which(colSums(NR) == 0))
  NR=NR[,!colnames(NR)%in%samples_exclude_int]
  indices_selected <- which(rowSums(NR, na.rm=TRUE) > 0)
  cat("Number of rows with depth in at least 1 sample: ", length(indices_selected)," out of ",nrow(NR)," total rows ( ",length(indices_selected)/nrow(NR),")\n")
  NR <- NR[indices_selected,]
  NV = fread(nv_path,data.table=F)
  rownames(NV)=NV[,1]
  NV=NV[,-1]
  samples_exclude_int <- names(which(colSums(NV) == 0))
  NV=NV[,!colnames(NV)%in%samples_exclude_int]
  NV <- NV[indices_selected,]
  samples=colnames(NV)
  Muts=rownames(NV)
}else{
  print("Please provide either NV and NR files or a path to CGPVaf output")
  break
}

# Guard: no variants have depth > 0 in any sample (e.g. empty chr)
if (length(Muts) == 0) {
  cat("No variants with depth > 0 — writing empty output and exiting.\n")
  empty_df <- data.frame(
    Chr = character(), Pos = character(), Ref = character(), Alt = character(),
    Mean_Depth = numeric(), Depth_filter = logical(),
    Germline_pval = numeric(), Germline_pval_log10 = numeric(),
    Rho = numeric(), Gender = character()
  )
  write.table(empty_df,
              paste0(output_dir, patient_ID, "_both_filtering_all.txt"),
              sep = "\t", append = FALSE, quote = FALSE,
              row.names = TRUE, col.names = TRUE)
  quit(save = "no", status = 0)
}

#Solution for when nr is lower thant total nv across
Tab_NR <- NR
tot_sum <- 0
for( i in 1:ncol(NV)){
  y <- NV[,i]
  x <- NR[,i]
  tot_sum <- tot_sum + length(which(x < y))
  x[which(x < y)] <- y[which(x < y)]
  Tab_NR[,i] <- x

}

NR <- Tab_NR

all_total <- nrow(NR)*ncol(NR)
cat("Number of positions in which NR is lower than NV: ",tot_sum," of total ",all_total," positions (",tot_sum/all_total,"}\n")


Muts_coord=matrix(ncol=4,unlist(strsplit(Muts,split="_")),byrow = T)
if(all(nchar(Muts_coord[,3])==1&nchar(Muts_coord[,4]))==1){
  mut_id="snv"
} else{
  if(all(nchar(Muts_coord[,3])>1|nchar(Muts_coord[,4])>1)){
    mut_id="indel"
  } else{
    mut_id="both"
  }
}
print(paste0("Mutations in data:", mut_id))

if(is.null(gender_passed)){

  XY_chromosomal = grepl("X|Y",Muts)
  autosomal = !XY_chromosomal
  xy_depth=mean(rowMeans(NR[XY_chromosomal,]))
  autosomal_depth=mean(rowMeans(NR[autosomal,]))

  gender='male'
  if(xy_depth>0.8*autosomal_depth) gender='female'
}else{

  gender = gender_passed

}

noCNVs=!samples%in%samples_with_CNVs

#----------------------------------
# Filtering
#----------------------------------
if(output_dir!="") system(paste0("mkdir -p ",output_dir))
print("Starting filtering...")

filter_df=as.data.frame(matrix(ncol=4,unlist(strsplit(rownames(NV),split="_")),byrow = T))
rownames(filter_df)=rownames(NV)
colnames(filter_df)=c("Chr","Pos","Ref","Alt")

filter_df$Mean_Depth=rowMeans(NR[,noCNVs])

if(gender=='male'){

  #Here determine the type of chr we are proccesing
  XY_chromosomal = grepl("X|Y",rownames(NR))
  if(length(which(XY_chromosomal)==TRUE)>0){

    filter_df$Depth_filter =  (rowMeans(NR[,noCNVs])>(min_cov/2)&rowMeans(NR[,noCNVs])<(max_cov/2))

  }else{

    filter_df$Depth_filter = (rowMeans(NR[,noCNVs])>min_cov&rowMeans(NR[,noCNVs])<max_cov)

  }

}else{
  filter_df$Depth_filter = rowMeans(NR)>min_cov&rowMeans(NR)<max_cov
}

# Filter out variants likely to be germline
germline_qval=exact.binomial(gender=gender,NV=NV[,noCNVs],NR=NR[,noCNVs],qval_return=T) 
filter_df$Germline_qval=germline_qval
q_temp <- p.adjust(filter_df$Germline_qval, method = "fdr")
filter_df$Germline=as.numeric(log10(q_temp)<log10(0.05))

print("Running beta-binomial on shared mutations...")
  

NR_flt=NR[filter_df$Germline&
            filter_df$Depth_filter,]
NV_flt=NV[filter_df$Germline&
            filter_df$Depth_filter,]


NR_flt_nonzero=NR_flt
NR_flt_nonzero[NR_flt_nonzero==0]=1

# Find shared variants and run beta-binomial filter  
shared_muts=rownames(NV_flt)[rowSums(NV_flt>0)>1]

cat("Running beta binomial on ",length(shared_muts)," mutations using ",ncores," cores\n")

filter_df$Beta_binomial=filter_df$Rho=NA

if(length(shared_muts)>0){
  if(ncores>1){
    rho_est=unlist(mclapply(shared_muts,function(x){
      estimateRho_gridml(NR_vec=as.numeric(NR_flt_nonzero[x,]),NV_vec=as.numeric(NV_flt[x,]))
    },mc.cores=ncores))
  }else{
    rho_est = beta.binom.filter(NR=NR_flt_nonzero[shared_muts,],NV=NV_flt[shared_muts,])
  }

  filter_df[shared_muts,"Rho"]=rho_est
  filter_df[shared_muts,"Beta_binomial"]=1

  if(mut_id=="snv")flt_rho=rho_est<snv_rho
  if(mut_id=="indel")flt_rho=rho_est<indel_rho
  if(mut_id=="both"){
    Muts_coord=matrix(ncol=4,unlist(strsplit(shared_muts,split="_")),byrow = T)
    is.indel=nchar(Muts_coord[,3])>1|nchar(Muts_coord[,4])>1
    flt_rho=(rho_est<indel_rho&is.indel)|(rho_est<snv_rho&!is.indel)
  }
  rho_filtered_out = shared_muts[flt_rho]
  filter_df[rho_filtered_out,"Beta_binomial"]=0
}

filter_df$Gender <- gender
filter_df$Germline_qval_log10 <- log10(filter_df$Germline_qval)
filter_df <- filter_df[,c("Chr","Pos","Ref","Alt","Mean_Depth","Depth_filter","Germline_qval","Germline_qval_log10","Rho","Gender")]
colnames(filter_df) <- c("Chr","Pos","Ref","Alt","Mean_Depth","Depth_filter","Germline_pval","Germline_pval_log10","Rho","Gender")
write.table(filter_df,paste0(output_dir,patient_ID,"_",mut_id,"_filtering_all.txt"),sep = "\t",append = FALSE,quote = FALSE,row.names = TRUE,col.names = TRUE)