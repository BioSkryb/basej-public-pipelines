
if(!require("optparse", character.only=T,quietly = T, warn.conflicts = F)){
  install.packages("optparse",repos = "http://cran.us.r-project.org")
  library("optparse", character.only=T,quietly = T, warn.conflicts = F)
}
#----------------------------------
# Input options
#----------------------------------
option_list = list(
  make_option(c("-i", "--donor_id"), action="store", default='Patient', type='character', help="Patient/donor ID to add to names of output files"),
  make_option(c("-v", "--input_nv"), action="store", default=NULL, type='character', help="Input NV matrix (rows are variants, columns are samples)"),
  make_option(c("-r", "--input_nr"), action="store", default=NULL, type='character', help="Input NR matrix (rows are variants, columns are samples)"),
  make_option(c("-c", "--cgpvaf_output"), action="store", default=NULL, type='character', help="CGPVaf output file, instead of NR/NV matrices - can be multiple files, i.e. indel and snv data for the same donor (comma-separated)"),
  make_option(c("-o", "--output_dir"), action="store", default="", type='character', help="Output directory for files"),
  make_option(c("-b", "--beta_binom_shared"), action="store", default=T, type='logical', help="Only run beta-binomial filter on shared mutations. If FALSE, run on all mutations, before germline/depth filtering"),
  make_option(c("-n", "--ncores"), action="store", default=1, type='numeric', help="Number of cores to use for the beta-binomial step"),
  make_option(c("--normal_flt"), action="store", default='PDv37is', type='character', help="Name of the dummy normal to exclude from cgpVAF output"),
  make_option(c("--snv_rho"), action="store", default=0.1, type='numeric', help="Rho value threshold for SNVs"),
  make_option(c("--indel_rho"), action="store", default=0.15, type='numeric', help="Rho value threshold for indels"),
  make_option(c("--min_cov"), action="store", default=10, type='numeric', help="Lower threshold for mean coverage across variant site"),
  make_option(c("--max_cov"), action="store", default=500, type='numeric', help="Upper threshold for mean coverage across variant site"),
  make_option(c("--only_snvs"), action="store", default=T, type='logical', help="If indel file is provided, only use SNVs to construct the tree (indels will still be mapped to branches)"),
  make_option(c("--split_trees"), action="store", default=T, type='logical', help="If both indels and SNVs are provided, plot trees separately for each."),
  make_option(c("--keep_ancestral"), action="store", default=F, type='logical', help="Keep an ancestral branch in the phylogeny for mutation mapping"),
  make_option(c("-x","--exclude_samples"), action="store", default=NULL, type='character', help="Option to manually exclude certain samples from the analysis, separate with a comma"),
  make_option(c("--cnv_samples"), action="store", default=NULL, type='character', help="Samples with CNVs, exclude from germline/depth-based filtering, separate with a comma"),
  make_option(c("--vaf_absent"), action="store", default=0.1, type='numeric', help="VAF threshold (autosomal) below which a variant is absent"),
  make_option(c("--vaf_present"), action="store", default=0.3, type='numeric', help="VAF threshold (autosomal) above which a variant is present"),
  make_option(c("--gender"), action="store", default=NULL, type='character', help="Override gender inference: 'male' or 'female'. If NULL, inferred from X/Y depth ratio."),
  make_option(c("-m", "--mixmodel"), action="store", default=F, type='logical', help="Use a binomial mixture model to filter out non-clonal samples?"),
  make_option(c("--min_clonal_mut"), action="store", default=35, type='numeric', help="If using binomial mixture model, minimum number of clonal mutations (in cluster higher than --VAF_treshold_mixmodel) needed to include sample."),
  make_option(c("-t", "--tree_mut_pval"), action="store", default=0.01, type='numeric', help="Pval threshold for treemut's mutation assignment"),
  make_option(c("-g", "--genotype_conv_prob"), action="store", default=F, type='logical', help="Use a binomial mixture model to filter out non-clonal samples?"),
  make_option(c("-p", "--min_pval_for_true_somatic"), action="store", default=0.05, type='numeric', help="Pval threshold for somatic presence if generating a probabilistic genotype matrix"),
  make_option(c("--min_variant_reads_shared"), action="store", default=2, type='numeric', help="Minimum variant reads used in generating a probabilistic genotype matrix"),
  make_option(c("--min_vaf_shared"), action="store", default=2, type='numeric', help="Minimum VAF used in generating a probabilistic genotype matrix"),
  make_option(c("--create_multi_tree"), action="store", default=T, type='logical', help="Convert dichotomous tree from MPBoot to polytomous tree"),
  make_option(c("--mpboot_path"), action="store", default="", type='character', help="Path to MPBoot executable"),
  make_option(c("--germline_cutoff"), action="store", default=-5, type='numeric', help="Log10 of germline qval cutoff"),
  make_option(c("--genomeFile"), action="store", default="/nfs/cancer_ref01/Homo_sapiens/37/genome.fa", type='character', help="Reference genome fasta for plotting mutational spectra"),
  make_option(c("--plot_spectra"), action="store", default=F, type='logical', help="Plot mutational spectra?"),
  make_option(c("--max_muts_plot"), action="store", default=5000, type='numeric', help="Maximum number of SNVs to plot in mutational spectra"),
  make_option(c("--lowVAF_filter"), action="store", default=0, type='numeric', help="Minimum VAF threshold to filter out subclonal variants. Disabled by default."),
  make_option(c("--lowVAF_filter_positive_samples"), action="store", default=0, type='numeric', help="Read number to apply exact binomial filter for samples with more than given number of reads. Disabled by default."),
  make_option(c("--VAF_treshold_mixmodel"), action="store", default=0.3, type='numeric', help="VAF threshold for the mixture modelling step to consider a sample clonal")
)
opt = parse_args(OptionParser(option_list=option_list, add_help_option=T))

print(opt)

dp_pos=opt$lowVAF_filter_positive_samples
ncores=opt$ncores
lowVAF_threshold=opt$lowVAF_filter
normal_flt=opt$normal_flt
snv_rho=opt$snv_rho
genomeFile=opt$genomeFile
plot_spectra=opt$plot_spectra
VAF_treshold=opt$VAF_treshold_mixmodel
indel_rho=opt$indel_rho
min_cov=opt$min_cov
max_cov=opt$max_cov
output_dir=opt$output_dir
only_snvs=opt$only_snvs
germline_cutoff=opt$germline_cutoff
if(is.null(opt$exclude_samples)) {samples_exclude=NULL} else {samples_exclude=unlist(strsplit(x=opt$exclude_samples,split = ","))}
if(is.null(opt$cnv_samples)) {samples_with_CNVs=NULL} else {samples_with_CNVs=unlist(strsplit(x=opt$cnv_samples,split = ","))}
if(is.null(opt$cgpvaf_output)) {cgpvaf_paths=NULL} else {cgpvaf_paths=unlist(strsplit(x=opt$cgpvaf_output,split = ","))}
keep_ancestral=opt$keep_ancestral
patient_ID=opt$donor_id
output_dir=opt$output_dir
nv_path=opt$input_nv
nr_path=opt$input_nr
max_muts_plot=opt$max_muts_plot
VAF_present=opt$vaf_present
VAF_absent=opt$vaf_absent
mixmodel=opt$mixmodel
split_trees=opt$split_trees
genotype_conv_prob=opt$genotype_conv_prob
min_pval_for_true_somatic_SHARED = opt$min_pval_for_true_somatic
min_variant_reads_SHARED=opt$min_variant_reads_shared
min_vaf_SHARED=opt$min_vaf_shared
tree_mut_pval=opt$tree_mut_pval
beta_binom_shared=opt$beta_binom_shared
create_multi_tree=opt$create_multi_tree
path_to_mpboot=opt$mpboot_path
min_clonal_mut=opt$min_clonal_mut

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

plot_spectrum = function(bed,save,add_to_title="",genomeFile = "/nfs/cancer_ref01/Homo_sapiens/37/genome.fa"){
  mutations=as.data.frame(bed)
  colnames(mutations) = c("chr","pos","ref","mut")
  mutations$pos=as.numeric(mutations$pos)
  mutations = mutations[(mutations$ref %in% c("A","C","G","T")) & (mutations$mut %in% c("A","C","G","T")) & mutations$chr %in% c(1:22,"X","Y"),]
  mutations$trinuc_ref = as.vector(scanFa(genomeFile, GRanges(mutations$chr, IRanges(as.numeric(mutations$pos)-1, 
                                                                                     as.numeric(mutations$pos)+1))))
  # 2. Annotating the mutation from the pyrimidine base
  ntcomp = c(T="A",G="C",C="G",A="T")
  mutations$sub = paste(mutations$ref,mutations$mut,sep=">")
  mutations$trinuc_ref_py = mutations$trinuc_ref
  for (j in 1:nrow(mutations)) {
    if (mutations$ref[j] %in% c("A","G")) { # Purine base
      mutations$sub[j] = paste(ntcomp[mutations$ref[j]],ntcomp[mutations$mut[j]],sep=">")
      mutations$trinuc_ref_py[j] = paste(ntcomp[rev(strsplit(mutations$trinuc_ref[j],split="")[[1]])],collapse="")
    }
  }
  
  # 3. Counting subs
  freqs = table(paste(mutations$sub,paste(substr(mutations$trinuc_ref_py,1,1),substr(mutations$trinuc_ref_py,3,3),sep="-"),sep=","))
  sub_vec = c("C>A","C>G","C>T","T>A","T>C","T>G")
  ctx_vec = paste(rep(c("A","C","G","T"),each=4),rep(c("A","C","G","T"),times=4),sep="-")
  full_vec = paste(rep(sub_vec,each=16),rep(ctx_vec,times=6),sep=",")
  freqs_full = freqs[full_vec]; freqs_full[is.na(freqs_full)] = 0; names(freqs_full) = full_vec
  
  xstr = paste(substr(full_vec,5,5), substr(full_vec,1,1), substr(full_vec,7,7), sep="")
  
  if(!is.null(save)) pdf(save,width=12,height=4)
  if(is.null(save)) dev.new(width=12,height=4)
  colvec = rep(c("dodgerblue","black","red","grey70","olivedrab3","plum2"),each=16)
  y = freqs_full; maxy = max(y)
  h = barplot(y, las=2, col=colvec, border=NA, ylim=c(0,maxy*1.5), space=1, cex.names=0.6, names.arg=xstr, ylab="Number mutations", main=paste0("Number of mutations: ",sum(freqs_full), add_to_title))
  for (j in 1:length(sub_vec)) {
    xpos = h[c((j-1)*16+1,j*16)]
    rect(xpos[1]-0.5, maxy*1.2, xpos[2]+0.5, maxy*1.3, border=NA, col=colvec[j*16])
    text(x=mean(xpos), y=maxy*1.3, pos=3, label=sub_vec[j])
  }    
  if(!is.null(save)) dev.off()
  
}

exact.binomial=function(gender,NV,NR,cutoff=-5,qval_return=F){
  # Function to filter out germline variants based on unmatched
  # variant calls of multiple samples from same individual (aggregate coverage
  # ideally >150 or so, but will work with less). NV is matrix of reads supporting
  # variants and NR the matrix with total depth (samples as columns, mutations rows,
  # with rownames as chr_pos_ref_alt or equivalent). Returns a logical vector,
  # TRUE if mutation is likely to be germline.

  XY_chromosomal = grepl("X|Y",rownames(NR))
  flag_is_xy <- FALSE
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

    if(flag_is_xy==FALSE){

      cat("Entered male autosomal binomial filtering. FLAG_XY: ",flag_is_xy,"\n")

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
  qval = p.adjust(pval,method="BH")
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

binom_pval_matrix = function(NV,NR,gender,qval_return=F) {
  NR_nonzero=NR
  NR_nonzero[NR_nonzero==0]=1
  pval_mat <- matrix(0, nrow = nrow(NV), ncol = ncol(NV))
  rownames(pval_mat)=rownames(NV)
  colnames(pval_mat)=colnames(NV)
  if(gender == "male") {
    for(i in 1:nrow(NV)) {
      for (j in 1:ncol(NV)) {
        if (!grepl("X|Y",rownames(NV)[1])) {pval_mat[i,j] <- binom.test(NV[i,j], NR_nonzero[i,j], p = 0.5, alternative = "less")$p.value}
        else {pval_mat[i,j] <- binom.test(NV[i,j], NR_nonzero[i,j], p = 0.95, alternative = "less")$p.value}
      }
    }
  } else if(gender == "female") {
    for(i in 1:nrow(NV)) {
      for (j in 1:ncol(NV)) {
        pval_mat[i,j] <- binom.test(NV[i,j], NR_nonzero[i,j], p = 0.5, alternative = "less")$p.value
      }
    }
  }
  if(qval_return){
    qval_mat=matrix(p.adjust(as.vector(pval_mat), method='BH'),ncol=ncol(pval_mat))
    rownames(qval_mat)=rownames(NV)
    colnames(qval_mat)=colnames(NV)
    return(qval_mat)
  }else{
    return(pval_mat)
  }
}


apply_mix_model=function(NV,NR,plot=T,min_clonal_mut_num=min_clonal_mut){
  peak_VAF=rep(0,ncol(NV))
  names(peak_VAF)=colnames(NV)
  autosomal=!grepl("X|Y",rownames(NV))
  for(s in colnames(NV)){
    muts_include=NV[,s]>3&autosomal
    if(sum(muts_include)>5){
      NV_vec=NV[muts_include,s]
      NR_vec=NR[muts_include,s]
      res=binom_mix(NV_vec,NR_vec,mode="Truncated",nrange=1:3)
      saveRDS(res,paste0(output_dir,s,"_binom_mix.Rdata"))
      
      if(plot){
        pdf(paste0(output_dir,s,"_binom_mix.pdf"))
        p=hist(NV_vec/NR_vec,breaks=20,xlim=c(0,1),col='gray',freq=F,xlab="Variant Allele Frequency",
               main=paste0(s,", (n=",length(NV_vec),")"))
        cols=c("red","blue","green","magenta","cyan")
        
        y_coord=max(p$density)-0.5
        y_intv=y_coord/5
        
        for (i in 1:res$n){
          depth=rpois(n=5000,lambda=median(NR_vec))
          sim_NV=unlist(lapply(depth,rbinom,n=1,prob=res$p[i]))
          sim_VAF=sim_NV/depth
          sim_VAF=sim_VAF[sim_NV>3]
          dens=density(sim_VAF)
          lines(x=dens$x,y=res$prop[i]*dens$y,lwd=2,lty='dashed',col=cols[i])
          y_coord=y_coord-y_intv/2
          text(y=y_coord,x=0.9,label=paste0("p1: ",round(res$p[i],digits=2)))
          segments(lwd=2,lty='dashed',col=cols[i],y0=y_coord+y_intv/4,x0=0.85,x1=0.95)
        }
        dev.off()
      }
      peak_VAF[s]=max(res$p[(res$prop*length(res$Which_cluster))>min_clonal_mut])
    }
  }
  return(peak_VAF)
}

add_ancestral_outgroup=function(tree,outgroup_name="Ancestral"){
  #This function adds the ancestral tip at the end
  tmp=tree$edge
  N=length(tree$tip.label)
  newroot=N+2
  renamedroot=N+3
  ancestral_tip=N+1
  tmp=ifelse(tmp>N,tmp+2,tmp)
  
  tree$edge=rbind(c(newroot,renamedroot),tmp,c(newroot,ancestral_tip))
  tree$edge.length=c(0,tree$edge.length,0)
  
  tree$tip.label=c(tree$tip.label,outgroup_name)
  tree$Nnode=tree$Nnode+1
  mode(tree$Nnode)="integer"
  mode(tree$edge)="integer"
  return(tree)
}

low_vaf_in_pos_samples = function(NR, NV, gender, define_pos = 3, qval_return = F) {
  pval=rep(0,nrow(NR))
  if(gender == "male") {
    for(n in 1:nrow(NR)) {
      NV_vec=NV[n,]
      NR_vec=NR[n,]
      if(any(NV_vec >= define_pos)){
        NV_vec_pos=NV_vec[which(NV_vec >= define_pos)]
        NR_vec_pos=NR_vec[which(NV_vec >= define_pos)]
        if (grepl("X|Y",rownames(NR)[n])) {
          pval[n]=binom.test(sum(NV_vec_pos), sum(NR_vec_pos), p = 0.95, alt = "less")$p.value
        } else {
          pval[n]=binom.test(sum(NV_vec_pos), sum(NR_vec_pos), p = 0.5, alt = "less")$p.value
        }
      }
    }
  } else if(gender == "female") {
    for(n in 1:nrow(NR)) {
      NV_vec=NV[n,]
      NR_vec=NR[n,]
      if(any(NV_vec >= define_pos)){
        NV_vec_pos=NV_vec[which(NV_vec >= define_pos)]
        NR_vec_pos=NR_vec[which(NV_vec >= define_pos)]
        pval[n]=binom.test(sum(NV_vec_pos), sum(NR_vec_pos), p = 0.5, alt = "less")$p.value
      }
    }
  }
  if(qval_return){
    return(p.adjust(pval,method="BH"))
  }else{
    return(pval)
  }
}

#----------------------------------
# Read in data
#----------------------------------
print("Reading in data...")

if(!is.null(cgpvaf_paths)){
  if(length(cgpvaf_paths)==1){
    data = fread(cgpvaf_paths,header=T,data.table=F)
    Muts = paste(data$Chrom,data$Pos,data$Ref,data$Alt,sep="_")
    NR = data[,grepl("DEP",colnames(data))&!grepl(paste(c(normal_flt,samples_exclude),collapse="|"),colnames(data))]
    NV = data[,grepl("MTR",colnames(data))&!grepl(paste(c(normal_flt,samples_exclude),collapse="|"),colnames(data))]
    rownames(NV)=rownames(NR)=Muts
    samples=colnames(NR)=colnames(NV)=gsub("_DEP","",colnames(NR))
  }else{
    NR=NV=Muts=c()
    for(n in 1:length(cgpvaf_paths)){
      data = fread(cgpvaf_paths[n],header=T,data.table=F)
      Muts = c(Muts,paste(data$Chrom,data$Pos,data$Ref,data$Alt,sep="_"))
      NR = rbind(NR,data[,grepl("DEP",colnames(data))&!grepl(paste(c(normal_flt,samples_exclude),collapse="|"),colnames(data))])
      NV = rbind(NV,data[,grepl("MTR",colnames(data))&!grepl(paste(c(normal_flt,samples_exclude),collapse="|"),colnames(data))])
    }
    rownames(NV)=rownames(NR)=Muts
    samples=colnames(NR)=colnames(NV)=gsub("_DEP","",colnames(NR))
  }
}else{    
  if(!is.null(nr_path)&!is.null(nv_path)){
    NR = fread(nr_path,data.table=F)
    rownames(NR)=NR[,1]
    NR=NR[,-1]
    samples_exclude <- names(which(colSums(NR) == 0))
    NR=NR[,!colnames(NR)%in%samples_exclude]
    NV = fread(nv_path,data.table=F)
    rownames(NV)=NV[,1]
    NV=NV[,-1]
    samples_exclude <- names(which(colSums(NV) == 0))
    NV=NV[,!colnames(NV)%in%samples_exclude]
    samples=colnames(NV)
    Muts=rownames(NV)
  }else{
    print("Please provide either NV and NR files or a path to CGPVaf output")
    break
  }
}

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

XY_chromosomal = grepl("X|Y",Muts)
autosomal = !XY_chromosomal

if(!is.null(opt$gender)){
  stopifnot(opt$gender %in% c("male","female"))
  gender=opt$gender
} else {
  if(any(XY_chromosomal) && any(autosomal)){
    xy_depth=mean(rowMeans(NR[XY_chromosomal,,drop=FALSE]))
    autosomal_depth=mean(rowMeans(NR[autosomal,,drop=FALSE]))
    gender='male'
    if(xy_depth>0.8*autosomal_depth) gender='female'
  } else {
    gender='female'
  }
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
# Filter out variant sites with high and low depth across samples
if(gender=='male'){
  if(any(XY_chromosomal)){
    filter_df$Depth_filter = (rowMeans(NR[,noCNVs])>(min_cov/2)&rowMeans(NR[,noCNVs])<(max_cov/2))
  }else{
    filter_df$Depth_filter = (rowMeans(NR[,noCNVs])>min_cov&rowMeans(NR[,noCNVs])<max_cov)
  }
}else{
  filter_df$Depth_filter = rowMeans(NR)>min_cov&rowMeans(NR)<max_cov
}

# Filter out variants likely to be germline
germline_qval=exact.binomial(gender=gender,NV=NV[,noCNVs],NR=NR[,noCNVs],qval_return=T) 
filter_df$Germline_qval=germline_qval
filter_df$Germline=as.numeric(log10(germline_qval)<germline_cutoff)


if(lowVAF_threshold>0){
  NR_nonzero=NR
  NR_nonzero[NR_nonzero==0]=1
  VAF=NV/NR_nonzero
  filter_df$lowVAF=rowSums(VAF>lowVAF_threshold)>0
}

if(beta_binom_shared){
  print("Running beta-binomial on shared mutations...")
  
  if(lowVAF_threshold>0){
    NR_flt=NR[filter_df$Germline&
                filter_df$Depth_filter&
                filter_df$lowVAF,]
    NV_flt=NV[filter_df$Germline&
                filter_df$Depth_filter&
                filter_df$lowVAF,]
  }else{
    NR_flt=NR[filter_df$Germline&
                filter_df$Depth_filter,]
    NV_flt=NV[filter_df$Germline&
                filter_df$Depth_filter,]
  }
  
  NR_flt_nonzero=NR_flt
  NR_flt_nonzero[NR_flt_nonzero==0]=1
  
  # Find shared variants and run beta-binomial filter  
  shared_muts=rownames(NV_flt)[rowSums(NV_flt>0)>1]
  
  if(ncores>1){
    rho_est=unlist(mclapply(shared_muts,function(x){
      estimateRho_gridml(NR_vec=as.numeric(NR_flt_nonzero[x,]),NV_vec=as.numeric(NV_flt[x,]))
    },mc.cores=ncores))
  }else{
    rho_est = beta.binom.filter(NR=NR_flt_nonzero[shared_muts,],NV=NV_flt[shared_muts,])
  }
  
  filter_df$Beta_binomial=filter_df$Rho=NA
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
  
  NR_filtered = NR_flt[!rownames(NR_flt)%in%rho_filtered_out,]
  NV_filtered = NV_flt[!rownames(NV_flt)%in%rho_filtered_out,]
  
  if(dp_pos>0){
    pval_vec=low_vaf_in_pos_samples(NR=NR_filtered,NV=NV_filtered,gender=gender,define_pos=dp_pos)
    filter_df$lowVAF_in_pos_samples=NA
    filter_df[rownames(NR_filtered),"lowVAF_in_pos_samples"]=pval_vec>1e-3
  }
  
}else{
  print("Running beta-binomial on ALL mutations...")
  
  if(ncores>1){
    rho_est=unlist(mclapply(1:nrow(NR),function(x){
      estimateRho_gridml(NR_vec=as.numeric(NR[x,]),NV_vec=as.numeric(NV[x,]))
    },mc.cores=ncores))
  }else{
    rho_est=beta.binom.filter(NR=NR, NV=NV)
  }
  
  filter_df$Rho=rho_est
  if(mut_id=="snv")filter_df$Beta_binomial=as.numeric(rho_est>snv_rho&!is.na(rho_est))
  if(mut_id=="indel")filter_df$Beta_binomial=as.numeric(rho_est>indel_rho&!is.na(rho_est))
  if(mut_id=="both"){
    is.indel=nchar(filter_df$Ref)>1|nchar(filter_df$Alt)>1
    filter_df$Beta_binomial=as.numeric(((rho_est>indel_rho&is.indel)|(rho_est>snv_rho&!is.indel))&!is.na(rho_est))
  }
  
  
  filter_names=c("Depth_filter","Germline","Beta_binomial")
  if(lowVAF_threshold>0){
    filter_names=c(filter_names,"lowVAF")
  }
  if(dp_pos>0){
    pval_vec=low_vaf_in_pos_samples(NR=NR,NV=NV,gender=gender,define_pos=dp_pos)
    filter_df$lowVAF_in_pos_samples_pval=pval
    filter_df$lowVAF_in_pos_samples=pval_vec>1e-3
    filter_names=c(filter_names,"lowVAF_in_pos_samples")
  }
  NV_filtered=NV[rowSums(filter_df[,filter_names])==length(filter_names),]
  NR_filtered=NR[rowSums(filter_df[,filter_names])==length(filter_names),]
}

write.table(NR_filtered,paste0(output_dir,patient_ID,"_",mut_id,"_NR_filtered_all.txt"))
write.table(NV_filtered,paste0(output_dir,patient_ID,"_",mut_id,"_NV_filtered_all.txt"))
write.table(filter_df,paste0(output_dir,patient_ID,"_",mut_id,"_filtering_all.txt"))
