# Include negative control into the analysis
# 3rd round: Biotin elution
# P6157 + P6159

# Set your WD
setwd('~/MP_RStudio_May2014/CellSurfaceExtractionKit/'); getwd();

# Read in the data
pos <- read.csv2("~/MP_RStudio_May2014/CellSurfaceExtractionKit/P6157/Isobar_P6157.csv", header=TRUE, stringsAsFactors=F)
neg <- read.csv2("~/MP_RStudio_May2014/CellSurfaceExtractionKit/P6159//Isobar_P6159.csv", header=TRUE, stringsAsFactors=F)
pos$protein.acs <-  gsub("[,-].*",'',pos$protein.acs)   # Labelled sample contains 729 proteins
neg$protein.acs <-  gsub("[,-].*",'',neg$protein.acs)   # Ctrl contains 74 proteins
dim(pos); dim(neg)

# Show the overlap graphically
library(gplots)
pdf("P6157_P6159_Overlap.pdf")
venn(list(pos=pos$protein.acs,
          neg=neg$protein.acs))
dev.off();

library(VennDiagram);
library(extrafont); 
grid.newpage();
venn.plot <- draw.pairwise.venn(729,74,63, c('Labelled sample','unlabelled control'), cex=rep(3.5,3),
                                col=c('plum1', 'plum3'), fill=c('plum1', 'plum3'), 
                                cat.pos=c(305,125), cat.dist=c(0.07,0.07), 
                                cat.fontfamily=rep('Franklin Gothic Book',2), cat.cex=rep(1.5,2), 
                                cat.col=c('plum1', 'plum3'),
                                scaled=TRUE, margin=c(0.08, 0.08), fontfamily=rep('Franklin Gothic Book',3));


# Create a combined dataframe
sc.df <- merge(pos[,c("protein.acs", "ID", "Description", "Gene", "global.spectra.count")], 
               neg[, c("protein.acs", "global.spectra.count")], 
               by="protein.acs")
names(sc.df) <- c("AC", "ID", "Description", "Gene", "GSC.pos", "GSC.neg")
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
sc.df <- na.omit(sc.df) # No NA values in this dataset
dim(sc.df)
write.csv2(sc.df, file="comb_P6157_P6159.csv")

# You can also start from here by reading in the data from the .csv file
sc.df <- read.csv2("~/MP_RStudio_May2014//CellSurfaceExtractionKit/comb_P6157_P6159.csv", header=T, stringsAsFactors=F)
dim(sc.df)

# Plot density
par(mfrow=c(1,1))
plot(density(sc.df$sqrt.ratio), xlab="Normalized spectral count ratios")


# =======================================================================
# Identify contaminants
cont <- sc.df[which(sc.df$ratio<=1), ]    # 7 contaminating proteins
clean.sf <- sc.df[-which(sc.df$ratio<=1), ]  # 56 proteins left
dim(sc.df); dim(cont); dim(clean.sf)

# Delete the contaminants from the sample dataset
length(intersect(cont$AC, pos$protein.acs))
pos.clean <- pos[grep(paste(cont$AC, collapse="|"), pos$protein.acs, invert=TRUE), ]; dim(pos.clean); dim(pos)
write.csv2(pos.clean, file="~/MP_RStudio_May2014/CellSurfaceExtractionKit/P6157/P6157_cleaned.csv")


# =======================================================================
# Show that within the proteins detected in lab sample and unlab ctrl, still PM proteins have a higher ratio pos/neg

# Version 1: use clean.sf (no proteins with ratio <= 1) - 56 proteins
# Read in the PM proteins identified in the sample
PM_6157 <- read.csv2("~/MP_RStudio_May2014/CellSurfaceExtractionKit/P6157/PM_6157.csv", header=TRUE, stringsAsFactors=F)
extr.PM <- clean.sf[grep(paste(PM_6157$protein.acs, collapse="|"), clean.sf$AC), ]
non.extr.PM <- clean.sf[grep(paste(PM_6157$protein.acs, collapse="|"), clean.sf$AC, invert=TRUE), ]
extr.PM <- na.omit(extr.PM)
non.extr.PM <- na.omit(non.extr.PM)
dim(extr.PM); dim(non.extr.PM)

pdf("~/MP_RStudio_May2014/CellSurfaceExtractionKit/pdf/Boxplot_P6157_P6159.scRatio.pdf")
boxplot(non.extr.PM$ratio, extr.PM$ratio, 
        names=c("non-PM", "PM"), 
        col=c("darkgrey", "red"), 
        ylab="Spectral Count Ratio")
dev.off();

# Version 2: on unfiltered data - 63 proteins
pdf("~/MP_RStudio_May2014/CellSurfaceExtractionKit/pdf/Boxplot_P6157_P6159.scRatioAllProteins.pdf")
boxplot(non.extr.PM.raw$ratio, extr.PM.raw$ratio, 
        names=c("non-PM", "PM"), 
        col=c("darkgrey", "red"), 
        ylab="Spectral Count Ratio")
dev.off();


# Identify "contaminants" in sample, i.e. proteins whose ratio is smaller than defined threshold
# Take different cutoffs and plot distribution
pdf("Distrib_scRatioDiffCutoffs_sqrt_P6157_P6159.pdf")
par(mfrow=c(2,3))
thres <- c(1,1.5,2,3,4,5)
for(thre in thres){
  cont <- sc.df[which(sc.df$ratio<=thre), ] 
  clean.sf <- sc.df[-which(sc.df$ratio<=thre), ]
  extr.PM <- clean.sf[grep(paste(PM_6157$protein.acs, collapse="|"), clean.sf$AC), ]
  non.extr.PM <- clean.sf[grep(paste(PM_6157$protein.acs, collapse="|"), clean.sf$AC, invert=TRUE), ]
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
extr.PM.raw <- sc.df[grep(paste(PM_6157$protein.acs, collapse="|"), sc.df$AC), ]
non.extr.PM.raw <- sc.df[grep(paste(PM_6157$protein.acs, collapse="|"), sc.df$AC, invert=TRUE), ]
extr.PM.raw <- na.omit(extr.PM.raw)
non.extr.PM.raw <- na.omit(non.extr.PM.raw)
dim(extr.PM.raw); dim(non.extr.PM.raw)

# log10 values of GSC contain -Inf values that need to be deleted
which(log10(extr.PM.raw$GSC.pos.norm)=="-Inf")
which(log10(extr.PM.raw$GSC.neg.norm)=="-Inf")
extr.PM.raw <- extr.PM.raw[-41,]
non.extr.PM.raw <- non.extr.PM.raw[-2,]

# Get an idea about data distribution to set axis limits for plotting

# Graphical representations
setwd("~/MP_RStudio_May2014/CellSurfaceExtractionKit/pdf")
pdf("sc_scatterplot_P6157_P6159.pdf")
par(mfrow=c(1,1))
plot(log10(extr.PM.raw$GSC.neg.norm), log10(extr.PM.raw$GSC.pos.norm), col="red", cex=.7, pch=21, bg="red",
     xlab="log10(sc(neg))", 
     ylab="log10(sc(pos))", 
     xlim=c(-0.8, 0.5), 
     ylim=c(-0.8,1.5),
     main="Cell Surface Protein Extraction Kit")
points(log10(non.extr.PM.raw$GSC.neg.norm), log10(non.extr.PM.raw$GSC.pos.norm), cex=.7, pch=21, bg="black")
abline(0,1,lty=2)
legend("bottomright", legend = c('PM', 'non-PM'), fill = c("red", "black"))
legend("topleft", legend=c("y=x"), lty=2)
dev.off()

# Selecting non-PM proteins with highest sc(pos)/sc(neg) ratio and label them in the scatterplot gets too messy
# Therefore sort the dataframe by ratio column
non.extr.PM.raw <- non.extr.PM.raw[order(non.extr.PM.raw[,10], decreasing=T), ]
setwd("~/MP_RStudio_May2014/CellSurfaceExtractionKit/")
write.csv2(non.extr.PM.raw, file="NonPM_Lab+Ctrl_P6157_P6159.csv")
