library(ctc)
library(DNAcopy) # Segmentation
library(inline)  # Use of c++
library(gplots)  # Visual plotting of tables
library(scales)
library(plyr)
library(ggplot2)
library(gridExtra)

args = commandArgs(trailingOnly=TRUE)
if (length(args)<8) {
	  stop("Seven arguments should be passed \n", call.=FALSE)
} 


dat <- args[1]
bm <- args[2]
gc_file <- args[3]
bounds_file <- args[4]
minPloidy <- as.numeric(args[5])
maxPloidy <- as.numeric(args[6])
minBinWidth <- as.numeric(args[7])
f <- as.numeric(args[8])
#stat is segmentation variable if it is passed then segmentation files was passed
#Here we create a dummy ref variable
#Default is 0
stat <- 0


#f denotes if facs file was provided
#default is 0

#bb denotes if bad bins should be clean
bb <- 0


### Load genome specific files ###
GC     = read.table(gc_file   , header=FALSE, sep="\t", as.is=TRUE)
loc    = read.table(bm                          , header=TRUE , sep="\t", as.is=TRUE)
bounds = read.table(bounds_file, header=FALSE, sep="\t")

# Load user data
raw    = read.table(dat, header=TRUE, sep="\t",check.names = FALSE)
ploidy = rbind(c(0,0), c(0,0))

#If segmentation was 1 or 2 then 
if (f == 1 | f == 2) {
  ploidy = read.table("ploidy.txt", header=FALSE, sep="\t", as.is=TRUE)  
}

#Remove bad bins
#Standardize this in function format
# Remove bad bins
if (bb)
{
  print("Removing bad bins...")
  badbins = read.table(paste(genome, "/badbins_", bm, sep=""), header=FALSE, sep="\t", as.is=TRUE)
  GC      = data.frame(GC[-badbins[,1], 1])
  loc     = loc[-badbins[,1], ]
  raw     = data.frame(raw[-badbins[,1], ])

  step  = 1
  chrom = loc[1,1]
  for (i in 1:nrow(loc))
  {
   if (loc[i,1] != chrom)
   {
     bounds[step,1] = chrom
     bounds[step,2] = i
     step           = step+1
     chrom          = loc[i,1]
    }
  }
}

# Initialize color palette
cp <- 3
colors     = matrix(0,3,2)
colors[1,] = c('goldenrod', 'darkmagenta')
colors[2,] = c('dodgerblue', 'darkorange')
#Row 3 in colors represent the deletion and duplication respectively
colors[3,] = c('#b2182b', '#2166ac')

#Initialize data structures
l            = dim(raw)[1] # Number of bins
w            = dim(raw)[2] # Number of cells
breaks       = matrix(0,l,w)
fixed        = matrix(0,l,w)
final        = matrix(0,l,w)
stats        = matrix(0,w,10)
pos          = cbind(c(1,bounds[,2]), c(bounds[,2], l))
# Initialize CN inference variables
CNgrid       = seq(minPloidy, maxPloidy, by=0.05)
n_ploidy     = length(CNgrid)  # Number of ploidy tests during CN inference
CNmult       = matrix(0,n_ploidy,w)
CNerror      = matrix(0,n_ploidy,w)
outerColsums = matrix(0,n_ploidy,w)

#Normalize cells
normal  = sweep(raw+1, 2, colMeans(raw+1), '/')
normal2 = normal
lab     = colnames(normal)

# Prepare statistics
rownames(stats) = lab
colnames(stats) = c("Reads", "Bins", "Mean", "Var", "Disp", "Min", "25th", "Median", "75th", "Max")

# Determine segmentation reference using dispersion (stat = 1) or reference sample (stat = 2)
if (stat == 1)
{
  F = normal[,which.min(apply(normal, 2, sd)/apply(normal,2,mean))[1]]
} else if (stat == 2) {
  R   = read.table(ref, header=TRUE, sep="\t", as.is=TRUE)
  low = lowess(GC[,1], log(R[,1]+0.001), f=0.05)
  app = approx(low$x, low$y, GC[,1])
  F   = exp(log(R[,1]) - app$y)
}

# Process all cells in the dataset
sink("results.txt")
cat(paste("Sample\tCopy_Number\tSoS_Predicted_Ploidy\tError_in_SoS_Approach\n", sep=""))
res_list_cnv <- list()
res_list_ploidy <- list()
for(k in 1:w){

   cat('===',k,'===\n')

  # Generate basic statistics
  stats[k,1]  = sum(raw[,k])
  stats[k,2]  = l
  stats[k,3]  = round(mean(raw[,k]), digits=2)
  stats[k,4]  = round(sd(raw[,k]), digits=2)
  stats[k,5]  = round(stats[k,4]/stats[k,3], digits=2)
  stats[k,6]  = min(raw[,k])
  stats[k,7]  = quantile(raw[,k], c(.25))[[1]]
  stats[k,8]  = median(raw[,k])
  stats[k,9]  = quantile(raw[,k], c(.75))[[1]]
  stats[k,10] = max(raw[,k])

  # ----------------------------------------------------------------------------
  # -- Segment data
  # ----------------------------------------------------------------------------

  # Calculate normalized for current cell (previous values of normal seem wrong)
  lowess.gc = function(jtkx, jtky) {
    jtklow = lowess(jtkx, log(jtky), f=0.05); 
    jtkz = approx(jtklow$x, jtklow$y, jtkx)
    return(exp(log(jtky) - jtkz$y))
  }
  normal[,k] = lowess.gc( GC[,1], (raw[,k]+1)/mean(raw[,k]+1) )

  # Compute log ratio between kth sample and reference
  if (stat == 0) {
    lr = log2(normal[,k])
  } else {
    lr = log2((normal[,k])/(F))
  }

  # Determine breakpoints and extract chrom/locations
  CNA.object   = CNA(genomdat = lr, chrom = loc[,1], maploc = as.numeric(loc[,2]), data.type = 'logratio')
  CNA.smoothed = smooth.CNA(CNA.object)
  segs         = segment(CNA.smoothed, verbose=0, min.width=minBinWidth)
  frag         = segs$output[,2:3]

  # Map breakpoints to kth sample
  len = dim(frag)[1]
  bps = array(0, len)
  for (j in 1:len)
    bps[j]=which((loc[,1]==frag[j,1]) & (as.numeric(loc[,2])==frag[j,2]))
  bps = sort(bps)
  bps[(len=len+1)] = l

  # Track global breakpoint locations
  breaks[bps,k] = 1

  # Modify bins to contain median read count/bin within each segment
  fixed[,k][1:bps[2]] = median(normal[,k][1:bps[2]])
  for(i in 2:(len-1))
    fixed[,k][bps[i]:(bps[i+1]-1)] = median(normal[,k][bps[i]:(bps[i+1]-1)])
  fixed[,k] = fixed[,k]/mean(fixed[,k])

  # ----------------------------------------------------------------------------
  # -- Determine Copy Number (SoS Method)
  # ----------------------------------------------------------------------------

  # Determine Copy Number     
  outerRaw         = fixed[,k] %o% CNgrid
  outerRound       = round(outerRaw)
  outerDiff        = (outerRaw - outerRound) ^ 2
  outerColsums[,k] = colSums(outerDiff, na.rm = FALSE, dims = 1)
  CNmult[,k]       = CNgrid[order(outerColsums[,k])]
  CNerror[,k]      = round(sort(outerColsums[,k]), digits=2)

  if (f == 0 | length(which(lab[k]==ploidy[,1]))==0 ) {
    CN = CNmult[1,k]
  } else if (f == 1) {
    CN = ploidy[which(lab[k]==ploidy[,1]),2]
    # If user specified FACS file, still calculate CNerror
    CNerror_facs = round( sort(colSums((round(fixed[,k] %o% c(CN)) - fixed[,k] %o% c(CN)) ^ 2, na.rm=FALSE, dims=1)), digits=2 )
  } else {
    estimate = ploidy[which(lab[k]==ploidy[,1]),2]
    CN = CNmult[which(abs(CNmult[,k] - estimate)<.4),k][1]
  }
  final[,k] = round(fixed[,k]*CN)

  #Create structures  for cnv plotting
  top=8
  rectangles1=data.frame(pos[seq(1,nrow(pos), 2),])
  rectangles2=data.frame(pos[seq(2,nrow(pos), 2),])
  clouds=data.frame(x=1:l, y=normal[,k]*CN)
  amp=data.frame(x=which(final[,k]>2), y=final[which(final[,k]>2),k])
  del=data.frame(x=which(final[,k]<2), y=final[which(final[,k]<2),k])
  flat=data.frame(x=which(final[,k]==2), y=final[which(final[,k]==2),k])
  anno=data.frame(x=(pos[,2]+pos[,1])/2, y=-top*.05, chrom=substring(c(as.character(bounds[,1]), "chrY"), 4 ,5))

  cnv_structure_plot <- list(top = top,
			     rectangles1 = rectangles1,
			     rectangles2 = rectangles2,
			     clouds = clouds,
			     amp = amp,
			     del = del,
			     flat = flat,
			     anno = anno)

  #Save structure in list
  res_list_cnv[[lab[k]]] <- cnv_structure_plot
  
  #Create list to hold estimate sample level CN estimate
  res_list_ploidy[[lab[k]]] <- CN

  #Perform within plotting here to verify 
  jpeg(filename=paste(lab[k], "_CN.jpeg", sep=""), width=3000, height=750)
    
  plot1 = ggplot() +
    geom_rect(data=rectangles1, aes(xmin=X1, xmax=X2, ymin=-top*.1, ymax=top), fill='gray85', alpha=0.75) +
    geom_rect(data=rectangles2, aes(xmin=X1, xmax=X2, ymin=-top*.1, ymax=top), fill='gray75', alpha=0.75) +
    geom_point(data=clouds, aes(x=x, y=y), color='gray45', size=3) +
    geom_point(data=flat, aes(x=x, y=y), size=4) +
    geom_point(data=amp, aes(x=x, y=y), size=4, color=colors[cp,1]) +
    geom_point(data=del, aes(x=x, y=y), size=4, color=colors[cp,2]) +
    geom_text(data=anno, aes(x=x, y=y, label=chrom), size=12) +
    scale_x_continuous(limits=c(0, l), expand = c(0, 0)) +
    scale_y_continuous(limits=c(-top*.1, top), expand = c(0, 0)) +
    labs(title=paste("Integer Copy Number Profile for Sample \"", lab[k], "\"\n Predicted Ploidy = ", CN, sep=""), x="Chromosome", y="Copy Number", size=16) +
    theme_minimal() +
    theme(plot.title=element_text(size=40, vjust=1.5)) +
    theme(axis.title.x=element_text(size=40, vjust=-.05), axis.title.y=element_text(size=40, vjust=.1)) +
    theme(axis.text=element_text(color="black", size=40), axis.ticks=element_line(color="black"))+
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.line.x = element_blank()) +
    theme(panel.background = element_rect(fill = 'gray90')) +
    theme(plot.margin=unit(c(.5,1,.5,1),"cm")) +
    theme(panel.grid.major.x = element_blank()) +
    geom_vline(xintercept = c(1, l), size=.5) +
    geom_hline(yintercept = c(-top*.1, top), size=.5)

    grid.arrange(plot1, ncol=1)
    
    dev.off()

  # Output results of CN calculations to file
  out=paste(lab[k], CN, paste(CNmult[,k], collapse= ","), paste(CNerror[,k], collapse= ","), sep="\t")
  cat(out, "\n")

}
sink()


# Store processed sample information
loc2=loc
loc2[,3]=loc2[,2]
pos = cbind(c(1,bounds[,2]), c(bounds[,2], l))

#
for (i in 1:nrow(pos))
{
  # If only 1 bin in a chromosome
  if( (pos[i,2] - pos[i,1]) == 0 ) {
    loc2[pos[i,1],1] = 1
  # If two bins.......
  } else if( (pos[i,2] - pos[i,1]) == 1 ) {
    loc2[pos[i,1],1] = 1
    loc2[pos[i,2],1] = loc2[pos[i,1],2] + 1
  } else {
    loc2[pos[i,1]:(pos[i,2]-1),2]=c(1,loc[pos[i,1]:(pos[i,2]-2),2]+1)
  }
}

# 
loc2[nrow(loc2),2]=loc2[nrow(loc2)-1,3]+1
colnames(loc2)=c("CHR","START", "END")

#Write final dirs
user_dir = "."
write.table(cbind(loc2,normal), file=paste(user_dir, "/SegNorm", sep=""), row.names=FALSE, col.names=c(colnames(loc2),lab), sep="\t", quote=FALSE)
write.table(cbind(loc2,fixed), file=paste(user_dir, "/SegFixed", sep=""), row.names=FALSE, col.names=c(colnames(loc2),lab), sep="\t", quote=FALSE)
write.table(cbind(loc2,final), file=paste(user_dir, "/SegCopy", sep=""), row.names=FALSE, col.names=c(colnames(loc2),lab), sep="\t", quote=FALSE)
write.table(cbind(loc2,breaks), file=paste(user_dir, "/SegBreaks", sep=""), row.names=FALSE, col.names=c(colnames(loc2),lab), sep="\t", quote=FALSE)
write.table(stats, file=paste(user_dir, "/SegStats", sep=""), sep="\t", quote=FALSE)

#Save structures
SegNorm <- cbind(loc2,normal)
SegFixed <- cbind(loc2,fixed)
SegCopy <- cbind(loc2,final)
SegBreaks <- cbind(loc2,breaks)

res_list <- list(SegNorm = SegNorm,
SegFixed = SegFixed,
SegCopy = SegCopy,
SegBreaks = SegBreaks,
CNV_structures = res_list_cnv,
CN = res_list_ploidy)

saveRDS(res_list ,file = "ginkgo_res.RDS")