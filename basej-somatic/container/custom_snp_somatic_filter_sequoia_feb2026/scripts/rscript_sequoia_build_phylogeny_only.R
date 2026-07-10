
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
  make_option(c("-o", "--output_dir"), action="store", default="", type='character', help="Output directory for files"),
  make_option(c("--only_snvs"), action="store", default=T, type='logical', help="If indel file is provided, only use SNVs to construct the tree (indels will still be mapped to branches)"),
  make_option(c("--split_trees"), action="store", default=T, type='logical', help="If both indels and SNVs are provided, plot trees separately for each."),
  make_option(c("--keep_ancestral"), action="store", default=F, type='logical', help="Keep an ancestral branch in the phylogeny for mutation mapping"),
  make_option(c("--vaf_absent"), action="store", default=0.1, type='numeric', help="VAF threshold (autosomal) below which a variant is absent"),
  make_option(c("--vaf_present"), action="store", default=0.3, type='numeric', help="VAF threshold (autosomal) above which a variant is present"),
  make_option(c("-t", "--tree_mut_pval"), action="store", default=0.01, type='numeric', help="Pval threshold for treemut's mutation assignment"),
  make_option(c("-g", "--genotype_conv_prob"), action="store", default=F, type='logical', help="Use a binomial mixture model to filter out non-clonal samples?"),
  make_option(c("-p", "--min_pval_for_true_somatic"), action="store", default=0.05, type='numeric', help="Pval threshold for somatic presence if generating a probabilistic genotype matrix"),
  make_option(c("--min_variant_reads_shared"), action="store", default=2, type='numeric', help="Minimum variant reads used in generating a probabilistic genotype matrix"),
  make_option(c("--min_vaf_shared"), action="store", default=2, type='numeric', help="Minimum VAF used in generating a probabilistic genotype matrix"),
  make_option(c("--create_multi_tree"), action="store", default=T, type='logical', help="Convert dichotomous tree from MPBoot to polytomous tree"),
  make_option(c("--gender"), action="store", default=NULL, type='character', help="Sample gender: 'male' or 'female' (required; used for chrX/Y VAF thresholds)"),
  make_option(c("--mpboot_path"), action="store", default="", type='character', help="Path to MPBoot executable"),
  make_option(c("--mpboot_bootstrap"), action="store", default=1000L, type="integer", help="MPBoot ultrafast bootstrap replicates (-bb); must be >= 1"),
  make_option(c("--input_genotype_bin"), action="store", default=NULL, type='character', help="Pre-computed genotype matrix TSV (rows=variants CHROM_POS_REF_ALT, cols=samples, values=0/0.5/1). When set, VAF discretization is skipped; NR and NV are still required for treemut branch assignment.")
)
opt = parse_args(OptionParser(option_list=option_list, add_help_option=T))

print(opt)

output_dir=opt$output_dir
only_snvs=opt$only_snvs
keep_ancestral=opt$keep_ancestral
patient_ID=opt$donor_id
nv_path=opt$input_nv
nr_path=opt$input_nr
VAF_present=opt$vaf_present
VAF_absent=opt$vaf_absent
split_trees=opt$split_trees
genotype_conv_prob=opt$genotype_conv_prob
min_pval_for_true_somatic_SHARED = opt$min_pval_for_true_somatic
min_variant_reads_SHARED=opt$min_variant_reads_shared
min_vaf_SHARED=opt$min_vaf_shared
tree_mut_pval=opt$tree_mut_pval
create_multi_tree=opt$create_multi_tree
path_to_mpboot=opt$mpboot_path
mpboot_bootstrap <- as.integer(opt$mpboot_bootstrap)
if (length(mpboot_bootstrap) != 1L || is.na(mpboot_bootstrap) || mpboot_bootstrap < 1L) {
  stop("--mpboot_bootstrap must be a positive integer, got: ", opt$mpboot_bootstrap)
}

#----------------------------------
# Load packages (install if they are not installed yet)
#----------------------------------
options(stringsAsFactors = F)
cran_packages=c("ape","seqinr","data.table")

for(package in cran_packages){
  if(!require(package, character.only=T,quietly = T, warn.conflicts = F)){
    install.packages(as.character(package),repos = "http://cran.us.r-project.org")
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

#----------------------------------
# Read in data
#----------------------------------
print("Reading in data...")

if(!is.null(nr_path)&!is.null(nv_path)){
  NR = fread(nr_path,data.table=F)
  rownames(NR)=NR[,1]
  NR=NR[,-1]
  NV = fread(nv_path,data.table=F)
  rownames(NV)=NV[,1]
  NV=NV[,-1]
  # Exclude samples with zero coverage in either matrix so NR and NV always
  # have identical columns (independent filtering could differ after indel
  # subsetting when only a handful of variants survive).
  samples_exclude <- union(names(which(colSums(NR) == 0)),
                           names(which(colSums(NV) == 0)))
  NR=NR[,!colnames(NR)%in%samples_exclude]
  NV=NV[,!colnames(NV)%in%samples_exclude]
  samples=colnames(NV)
  Muts=rownames(NV)
}else{
  print("Please provide both NV and NR files via --input_nv and --input_nr")
  break
}

# The input matrices are the already-filtered NR/NV
NV_filtered=NV
NR_filtered=NR

if (is.null(opt$gender)) stop("--gender is required: specify 'male' or 'female'")
gender <- tolower(opt$gender)
if (!gender %in% c("male", "female")) stop(sprintf("Invalid --gender value '%s': must be 'male' or 'female'", opt$gender))

if(output_dir!="") system(paste0("mkdir -p ",output_dir))

#----------------------------------
# Genotype matrix: load pre-built or discretize from NV/NR
#----------------------------------
if (!is.null(opt$input_genotype_bin)) {
  # Pre-built genotype_bin mode.
  # NR/NV are still used by treemut for branch assignment; only the VAF
  # discretization step is skipped.  Values in the file must be 0, 0.5, or 1.
  print("Loading pre-built genotype matrix (skipping VAF discretization)...")
  gb_raw <- fread(opt$input_genotype_bin, data.table = FALSE)
  rownames(gb_raw) <- gb_raw[, 1]
  genotype_bin <- as.matrix(gb_raw[, -1, drop = FALSE])
  mode(genotype_bin) <- "numeric"

  common_samples <- intersect(colnames(genotype_bin), samples)
  if (length(common_samples) == 0) {
    message("--input_genotype_bin sample names (first 10): ",
            paste(head(colnames(genotype_bin), 10), collapse = ", "))
    message("NR/NV matrix sample names (first 10): ",
            paste(head(samples, 10), collapse = ", "))
    stop("--input_genotype_bin: no overlapping samples with NR/NV matrices. ",
         "Check that column names in the genotype_bin TSV match the sample names ",
         "used in the NR/NV matrices (derived from VCF filenames with ",
         "_somatic_annotated[_filtered].vcf.gz stripped).")
  }
  common_muts <- intersect(rownames(genotype_bin), rownames(NR_filtered))
  if (length(common_muts) == 0) {
    message("--input_genotype_bin variant IDs (first 5): ",
            paste(head(rownames(genotype_bin), 5), collapse = ", "))
    message("NR/NV matrix variant IDs (first 5): ",
            paste(head(rownames(NR_filtered), 5), collapse = ", "))
    stop("--input_genotype_bin: no overlapping variants with NR/NV matrices. ",
         "Variant IDs must be CHROM_POS_REF_ALT (e.g. chr1_12345_A_T).")
  }

  genotype_bin <- genotype_bin[common_muts, common_samples, drop = FALSE]
  NR_filtered  <- NR_filtered[common_muts, common_samples, drop = FALSE]
  NV_filtered  <- NV_filtered[common_muts, common_samples, drop = FALSE]
  samples      <- common_samples
  print(paste0("Pre-built genotype_bin: ", nrow(genotype_bin), " variants x ", ncol(genotype_bin), " samples after intersection"))
} else {
  # Standard mode: discretize VAF into 0 / 0.5 / 1
  print("Constructing a fasta file...")
  NR_flt_nonzero=NR_filtered
  NR_flt_nonzero[NR_flt_nonzero==0]=1
  XY_chromosomal=grepl("X|Y",rownames(NR_filtered))
  autosomal=!XY_chromosomal

  if(genotype_conv_prob){
    pval_matrix=binom_pval_matrix(NR=NR_filtered, NV=NV_filtered,gender=gender)
    if(!is.na(min_variant_reads_SHARED)) {min_variant_reads_mat <- NV_filtered >= min_variant_reads_SHARED} else {min_variant_reads_mat=1}
    if(!is.na(min_pval_for_true_somatic_SHARED)) {min_pval_for_true_somatic_mat <- pval_matrix > min_pval_for_true_somatic_SHARED} else {min_pval_for_true_somatic_mat=1}
    if(!is.na(min_vaf_SHARED[1]) & gender=="female") {
      min_vaf_mat <- NV_filtered/NR_flt_nonzero>min_vaf_SHARED[1]
    } else if(!is.na(min_vaf_SHARED) & gender=="male") {
      min_vaf_mat=matrix(0,ncol=ncol(NV_filtered),nrow=nrow(NV_filtered))
      min_vaf_mat[XY_chromosomal,]=NV_filtered[XY_chromosomal,]/NR_flt_nonzero[XY_chromosomal,] > min_vaf_SHARED[2]
      min_vaf_mat[autosomal,]=NV_filtered[autosomal,]/NR_flt_nonzero[!autosomal,] > min_vaf_SHARED[1]
    } else {min_vaf_mat=1}
    genotype_bin = min_variant_reads_mat * min_pval_for_true_somatic_mat * min_vaf_mat
    #Select the "not sure" samples by setting genotype to 0.5.  THIS IS THE ONLY SLIGHTLY OPAQUE BIT OF THIS FUNCTION - SET EMPIRICALLY FROM EXPERIMENTATION.
    genotype_bin[NV_filtered > 0 & pval_matrix > 0.01 & genotype_bin != 1] <- 0.5 #If have any mutant reads, set as "?" as long as p-value > 0.01
    genotype_bin[NV_filtered >= 3 & pval_matrix > 0.001 & genotype_bin != 1] <- 0.5 #If have high numbers of mutant reads, should set as "?" even if incompatible p-value (may be biased sequencing)
    genotype_bin[(NV_filtered == 0) & (pval_matrix > 0.05)] <- 0.5 #Essentially if inadequate depth to exclude mutation, even if no variant reads
    # mut_id needed for pval matrix filename — derive here before the write
    Muts_coord_tmp=matrix(ncol=4,unlist(strsplit(rownames(NR_filtered),split="_")),byrow=T)
    mut_id_tmp=if(all(nchar(Muts_coord_tmp[,3])==1&nchar(Muts_coord_tmp[,4])==1)) "snv" else if(all(nchar(Muts_coord_tmp[,3])>1|nchar(Muts_coord_tmp[,4])>1)) "indel" else "both"
    write.table(pval_matrix,paste0(output_dir,patient_ID,"_",mut_id_tmp,"_filtered_binom_pval_mat.txt"))
  }else{
    genotype_bin=as.matrix(NV_filtered/NR_flt_nonzero)
    if(gender=="male"){
      genotype_bin[autosomal,][genotype_bin[autosomal,]<VAF_absent]=0
      genotype_bin[autosomal,][genotype_bin[autosomal,]>=VAF_present]=1
      genotype_bin[XY_chromosomal,][genotype_bin[XY_chromosomal,]<(2*VAF_absent)]=0
      genotype_bin[XY_chromosomal,][genotype_bin[XY_chromosomal,]>=(2*VAF_present)]=1
      genotype_bin[genotype_bin>0&genotype_bin<1]=0.5
    }
    if(gender=="female"){
      genotype_bin[genotype_bin<VAF_absent]=0
      genotype_bin[genotype_bin>=VAF_present]=1
      genotype_bin[genotype_bin>0&genotype_bin<1]=0.5
    }
  }
}

# mut_id — always derived from genotype_bin rows (consistent in both modes)
Muts_coord=matrix(ncol=4,unlist(strsplit(rownames(genotype_bin),split="_")),byrow = T)
if(all(nchar(Muts_coord[,3])==1&nchar(Muts_coord[,4])==1)){
  mut_id="snv"
} else{
  if(all(nchar(Muts_coord[,3])>1|nchar(Muts_coord[,4])>1)){
    mut_id="indel"
  } else{
    mut_id="both"
  }
}
print(paste0("Mutations in data:", mut_id))

present_vars_full=rowSums(genotype_bin>0)>0

if(only_snvs){
  Muts_coord=matrix(ncol=4,unlist(strsplit(rownames(genotype_bin),split="_")),byrow = T)
  is.indel=nchar(Muts_coord[,3])>1|nchar(Muts_coord[,4])>1
  genotype_bin=genotype_bin[!is.indel,]
}

#Create dummy fasta consisting of As (WT) and Ts (Mutant)
Ref = rep("A",nrow(genotype_bin))
Alt = rep("T",nrow(genotype_bin))
dna_strings = list()
dna_strings[1]=paste(Ref,sep="",collapse="") #Ancestral sample
for (n in 1:ncol(genotype_bin)){
  Mutations = Ref
  Mutations[genotype_bin[,n]==0.5] = '?'
  Mutations[genotype_bin[,n]==1] = Alt[genotype_bin[,n]==1]
  dna_string = paste(Mutations,sep="",collapse="")
  dna_strings[n+1]=dna_string
}

names(dna_strings)=c("Ancestral",colnames(genotype_bin))
# Short SC* tip labels for MPBoot FASTA (avoids long/special sample IDs in the alignment)
original_tip_names <- names(dna_strings)
mpboot_tip_rename_map <- data.frame(
  OriginalName = original_tip_names,
  NewName = paste0("SC", seq_along(original_tip_names)),
  stringsAsFactors = FALSE
)
write.table(
  mpboot_tip_rename_map,
  file = paste0(output_dir, patient_ID, "_", mut_id, "_mpboot_tip_rename_map.tsv"),
  quote = FALSE, row.names = FALSE, sep = "\t"
)
names(dna_strings) <- mpboot_tip_rename_map$NewName
require(seqinr)
write.fasta(dna_strings, names=names(dna_strings),paste0(output_dir,patient_ID,"_",mut_id,"_for_MPBoot.fa"))

#----------------------------------
# Build tree with MPBoot
#----------------------------------
print("Building a tree...")

system(paste0(path_to_mpboot,"mpboot-avx -s ",output_dir,patient_ID,"_",mut_id,"_for_MPBoot.fa -bb ", mpboot_bootstrap),ignore.stdout = T)

#----------------------------------
# Map Mutations on Tree using treemut
#----------------------------------

print("Mapping mutations...")
print("Assigning mutation without an ancestral branch")

tree=read.tree(paste0(output_dir,patient_ID,"_",mut_id,"_for_MPBoot.fa.treefile"))
# Restore biosample-style tip labels (FASTA used SC1, SC2, … for MPBoot)
tip_idx <- match(tree$tip.label, mpboot_tip_rename_map$NewName)
if (anyNA(tip_idx)) {
  stop("MPBoot treefile has tip labels not found in mpboot_tip_rename_map: ",
       paste(tree$tip.label[is.na(tip_idx)], collapse = ", "))
}
tree$tip.label <- mpboot_tip_rename_map$OriginalName[tip_idx]
tree=drop.tip(tree,"Ancestral")
if(!keep_ancestral){
  tree$edge.length=rep(1,nrow(tree$edge))
  NR_tree=NR_filtered[present_vars_full,]
  NV_tree=NV_filtered[present_vars_full,]
  res=assign_to_tree(tree,
                     mtr=as.matrix(NV_tree),
                     dep=as.matrix(NR_tree))
}else{
  tree <- add_ancestral_outgroup(tree) #Re add the ancestral outgroup after making tree dichotomous - avoids the random way that baseline polytomy is resolved
  tree$edge.length = rep(1, nrow(tree$edge))

  NR_tree=NR_filtered[present_vars_full,]
  NR_tree$Ancestral=30
  NV_tree=NV_filtered[present_vars_full,]
  NV_tree$Ancestral=0

  p.error = rep(0.01,ncol(NV_tree))
  p.error[colnames(NV_tree)=="Ancestral"]=1e-6
  res=assign_to_tree(tree,
                     mtr=as.matrix(NV_tree),
                     dep=as.matrix(NR_tree),
                     error_rate = p.error)
}

edge_length_nonzero = table(res$summary$edge_ml[res$summary$p_else_where<tree_mut_pval])
edge_length = rep(0,nrow(tree$edge))
names(edge_length)=1:nrow(tree$edge)
edge_length[names(edge_length_nonzero)]=edge_length_nonzero
tree$edge.length=as.numeric(edge_length)

if(create_multi_tree){
  print("Converting to a multi-furcating tree structure")
  if(keep_ancestral) {
    #Maintain the dichotomy with the ancestral branch
    ROOT=tree$edge[1,1]
    current_length<-tree$edge.length[tree$edge[,1]==ROOT & tree$edge[,2]!=which(tree$tip.label=="Ancestral")]
    new_length<-ifelse(current_length==0,1,current_length)
    tree$edge.length[tree$edge[,1]==ROOT & tree$edge[,2]!=which(tree$tip.label=="Ancestral")]<-new_length
  }
  tree<-di2multi(tree) #Now make tree multifurcating
  #Re-run the mutation assignment algorithm from the new tree
  res=assign_to_tree(tree,
                     mtr=as.matrix(NV_tree),
                     dep=as.matrix(NR_tree))
  edge_length_nonzero = table(res$summary$edge_ml[res$summary$p_else_where<tree_mut_pval])
  edge_length = rep(0,nrow(tree$edge))
  names(edge_length)=1:nrow(tree$edge)
  edge_length[names(edge_length_nonzero)]=edge_length_nonzero
  tree$edge.length=as.numeric(edge_length)
}

saveRDS(res,paste0(output_dir,patient_ID,"_",mut_id,"_assigned_to_tree.Rdata"))
write.tree(tree, paste0(output_dir,patient_ID,"_",mut_id,"_tree_with_branch_length_selectedscheme.tree"))

if(split_trees&mut_id=="both"){
  Muts_coord=matrix(ncol=4,unlist(strsplit(rownames(NV_filtered)[present_vars_full],split="_")),byrow = T)
  is.indel=nchar(Muts_coord[,3])>1|nchar(Muts_coord[,4])>1

  edge_length_nonzero = table(res$summary$edge_ml[res$summary$p_else_where<tree_mut_pval&!is.indel])
  edge_length = rep(0,nrow(tree$edge))
  names(edge_length)=1:nrow(tree$edge)
  edge_length[names(edge_length_nonzero)]=edge_length_nonzero
  tree$edge.length=as.numeric(edge_length)
  pdf(paste0(output_dir,patient_ID,"_snv_tree_with_branch_length.pdf"))
  plot(tree)
  axisPhylo(side = 1,backward=F)
  dev.off()
  write.tree(tree, paste0(output_dir,patient_ID,"_snv_tree_with_branch_length_selectedscheme.tree"))

  edge_length_nonzero = table(res$summary$edge_ml[res$summary$p_else_where<tree_mut_pval&is.indel])
  edge_length = rep(0,nrow(tree$edge))
  names(edge_length)=1:nrow(tree$edge)
  edge_length[names(edge_length_nonzero)]=edge_length_nonzero
  tree$edge.length=as.numeric(edge_length)
  pdf(paste0(output_dir,patient_ID,"_indel_tree_with_branch_length.pdf"))
  plot(tree)
  axisPhylo(side = 1,backward=F)
  dev.off()
  write.tree(tree, paste0(output_dir,patient_ID,"_indel_tree_with_branch_length_selectedscheme.tree"))

}else{
  pdf(paste0(output_dir,patient_ID,"_",mut_id,"_tree_with_branch_length.pdf"))
  plot(tree)
  axisPhylo(side = 1,backward=F)
  dev.off()

  tree_collapsed=tree
  tree_collapsed$edge.length=rep(1,nrow(tree_collapsed$edge))
  pdf(paste0(output_dir,patient_ID,"_",mut_id,"_tree_with_equal_branch_length.pdf"))
  plot(tree_collapsed)
  dev.off()
}

Mutations_per_branch=as.data.frame(matrix(ncol=4,unlist(strsplit(rownames(NR_tree),split="_")),byrow = T))
colnames(Mutations_per_branch)=c("Chr","Pos","Ref","Alt")
Mutations_per_branch$Branch = tree$edge[res$summary$edge_ml,2]
Mutations_per_branch=Mutations_per_branch[res$summary$p_else_where<tree_mut_pval,]
Mutations_per_branch$Patient = patient_ID
Mutations_per_branch$SampleID = paste(patient_ID,Mutations_per_branch$Branch,sep="_")
write.table(Mutations_per_branch,paste0(output_dir,patient_ID,"_",mut_id,"_assigned_to_branches.txt"),quote=F,row.names=F,sep="\t")

# Full placement table (one row per variant in NR_tree/NV_tree passed to assign_to_tree), including p_else_where
if (nrow(res$summary) != nrow(NR_tree)) {
  stop(sprintf("Placement summary length mismatch: res$summary has %d rows, NR_tree has %d",
               nrow(res$summary), nrow(NR_tree)))
}
placed_variants_all <- as.data.frame(matrix(ncol = 4,
  unlist(strsplit(rownames(NR_tree), split = "_")), byrow = TRUE), stringsAsFactors = FALSE)
colnames(placed_variants_all) <- c("Chr", "Pos", "Ref", "Alt")
placed_variants_all$edge_ml <- res$summary$edge_ml
placed_variants_all$Branch <- tree$edge[res$summary$edge_ml, 2]
placed_variants_all$p_else_where <- res$summary$p_else_where
placed_variants_all$pass_tree_mut_pval <- placed_variants_all$p_else_where < tree_mut_pval
placed_variants_all$Patient <- patient_ID
placed_variants_all$SampleID <- paste(patient_ID, placed_variants_all$Branch, sep = "_")
write.table(placed_variants_all,
            paste0(output_dir, patient_ID, "_", mut_id, "_placed_variants_all.tsv"),
            quote = FALSE, row.names = FALSE, sep = "\t")
saveRDS(placed_variants_all,
        paste0(output_dir, patient_ID, "_", mut_id, "_placed_variants_all.rds"))
