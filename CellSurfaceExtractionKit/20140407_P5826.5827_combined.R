# Include negative control into the analysis

# Set your WD
setwd('~/MP_RStudio_May2014/CellSurfaceExtractionKit/'); getwd();

# Read in the data
pos <- read.csv2("~/MP_RStudio_May2014/CellSurfaceExtractionKit/P5826/20140325_Isobar_P5826.csv", header=TRUE, stringsAsFactors=F)
neg <- read.csv2("~/MP_RStudio_May2014/CellSurfaceExtractionKit/P5827/20140325_Isobar_P5827.csv", header=TRUE, stringsAsFactors=F)
pos$protein.acs <-  gsub("[,-].*",'',pos$protein.acs)   # Labelled sample contains 966 proteins
neg$protein.acs <-  gsub("[,-].*",'',neg$protein.acs)   # Ctrl contains 606 proteins
dim(pos); dim(neg)

# Show the overlap graphically
library(gplots)
pdf("P5826_P5827_Overlap.pdf")
venn(list(pos=pos$protein.acs,
          neg=neg$protein.acs))
dev.off();

# Create a combined dataframe
sc.df <- merge(pos[,c("protein.acs", "ID", "Description", "Gene", "global.spectra.count")], 
               neg[, c("ID", "global.spectra.count")], 
               by="ID")
names(sc.df) <- c('ID', "AC", "Description", "Gene", "GSC.pos", "GSC.neg")
sc.df[["AC"]] <- gsub("[,-].*",'',sc.df[['AC']])    # Combined dataframe contains 499 proteins
dim(sc.df)

# Search for proteins that are detected with same GSC in positiv, labeled sample and negative ctrl
stable <- sc.df[which(sc.df$GSC.pos==sc.df$GSC.neg), ]
sum.stable.pos <- sum(stable$GSC.pos)
sum.stable.neg <- sum(stable$GSC.neg)

# Normalize the data to the sum of GSCs for stable subset
GSC.pos.norm <- sc.df$GSC.pos/sum.stable.pos
GSC.neg.norm <- sc.df$GSC.neg/sum.stable.neg
sc.df <- cbind(sc.df, GSC.pos.norm, GSC.neg.norm)

# Calculate the ratio between sample and neg. ctrl
ratio <- sc.df$GSC.pos.norm/sc.df$GSC.neg.norm
sqrt.ratio <- sqrt(ratio)
sc.df <- cbind(sc.df, ratio, sqrt.ratio)

# Clean the data frame of NA values
sc.df <- na.omit(sc.df) # 493 proteins left
dim(sc.df)
write.csv2(sc.df, file="comb_P5826_P5827.csv")

# You can also start from here by reading in the data from the .csv file
sc.df <- read.csv2("/slow/home/khoermann/MP/CellSurfaceExtractionKit/comb_P5826_P5827.csv", header=T, stringsAsFactors=F)
dim(sc.df)

# Plot density
par(mfrow=c(1,1))
plot(density(sc.df$sqrt.ratio), xlab="Normalized spectral count ratios")


# =======================================================================
# Identify contaminants
cont <- sc.df[which(sc.df$ratio<=1), ]    # 193 contaminating proteins
clean.sf <- sc.df[-which(sc.df$ratio<=1), ]  # 300 proteins left
dim(sc.df); dim(cont); dim(clean.sf)

# Delete the contaminants from the sample dataset
length(intersect(cont$AC, pos$protein.acs))
pos.clean <- pos[grep(paste(cont$AC, collapse="|"), pos$protein.acs, invert=TRUE), ]; dim(pos.clean); dim(pos)
write.csv2(pos.clean, file="~/MP_RStudio_May2014/CellSurfaceExtractionKit/P5826/P5826_cleaned.csv")

# =======================================================================
# Read in the PM proteins identified in the sample
PM_5826 <- read.csv2("~/MP_RStudio_May2014/CellSurfaceExtractionKit/P5826/PM_5826.csv", header=TRUE, stringsAsFactors=F)
extr.PM <- clean.sf[grep(paste(PM_5826$protein.acs, collapse="|"), clean.sf$AC), ]
non.extr.PM <- clean.sf[grep(paste(PM_5826$protein.acs, collapse="|"), clean.sf$AC, invert=TRUE), ]
extr.PM <- na.omit(extr.PM)
non.extr.PM <- na.omit(non.extr.PM)
dim(extr.PM); dim(non.extr.PM)

boxplot(non.extr.PM$ratio, extr.PM$ratio, names=c("non PM", "PM"))

# Identify "contaminants" in sample, i.e. proteins whose ratio is smaller than defined threshold
# Take different cutoffs and plot distribution
pdf("Distrib_scRatioDiffCutoffs_sqrt_P5826_P5827.pdf")
par(mfrow=c(2,3))
thres <- c(1,1.5,2,3,4,5)
for(thre in thres){
  cont <- sc.df[which(sc.df$ratio<=thre), ] 
  clean.sf <- sc.df[-which(sc.df$ratio<=thre), ]
  extr.PM <- clean.sf[grep(paste(PM_5826$protein.acs, collapse="|"), clean.sf$AC), ]
  non.extr.PM <- clean.sf[grep(paste(PM_5826$protein.acs, collapse="|"), clean.sf$AC, invert=TRUE), ]
  extr.PM <- na.omit(extr.PM)
  non.extr.PM <- na.omit(non.extr.PM)
  dens.PM <- density(extr.PM$sqrt.ratio)
  plot(density(non.extr.PM$sqrt.ratio), 
       xlim=c(0,7.5), ylim=c(0,1.5),
       xlab="Square root of normalized spectral count ratios", 
       main=paste("ratio cutoff:", thre, sep=" "))
  lines(dens.PM, col='red')
  legend("topright", legend = c('PM', 'non-PM'), fill = c("red", "black"))
}
dev.off()

# ===============================================================================
# For the scatterplot, it's better to leave the contaminants in
extr.PM.raw <- sc.df[grep(paste(PM_5826$protein.acs, collapse="|"), sc.df$AC), ]
non.extr.PM.raw <- sc.df[grep(paste(PM_5826$protein.acs, collapse="|"), sc.df$AC, invert=TRUE), ]
extr.PM.raw <- na.omit(extr.PM.raw)
non.extr.PM.raw <- na.omit(non.extr.PM.raw)
dim(extr.PM.raw); dim(non.extr.PM.raw)

# Graphical representations
setwd("~/MP_RStudio_May2014/CellSurfaceExtractionKit/pdf")
pdf("sc_scatterplot_P5826_P5827.pdf")
par(mfrow=c(1,1))
plot(log10(extr.PM.raw$GSC.neg.norm), log10(extr.PM.raw$GSC.pos.norm), col="red", cex=.7,
     xlab="log10(sc(neg))", 
     ylab="log10(sc(pos))", 
     main="Cell Surface Protein Extraction Kit")
points(log10(non.extr.PM.raw$GSC.neg.norm), log10(non.extr.PM.raw$GSC.pos.norm), cex=.7)
abline(0,1,lty=2)
legend("bottomright", legend = c('PM', 'non-PM'), fill = c("red", "black"))
legend("topleft", legend=c("y=x"), lty=2)
dev.off()

# Selcting non-PM proteins with highest sc(pos)/sc(neg) ratio and label them in the scatterplot gets too messy
# Therefore sort the dataframe by ratio column
non.extr.PM.raw <- non.extr.PM.raw[order(non.extr.PM.raw[,10], decreasing=T), ]
# Identify those where the ratio is bigger than 2
ind <- which(non.extr.PM.raw$ratio>2)
extreme.ratio.non.PM <- non.extr.PM.raw[ind, ]
setwd("~/MP_RStudio_May2014/CellSurfaceExtractionKit/")
write.xlsx2(extreme.ratio.non.PM, file="Extreme_Ratios_nonPM_cutoff2.xlsx")
