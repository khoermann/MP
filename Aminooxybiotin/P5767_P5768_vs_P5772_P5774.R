# Include negative control into the analysis

# Set your WD
setwd('~/MP_RStudio_May2014/Aminooxybiotin/'); getwd();

# Read in the data
pos <- read.csv2("~/MP_RStudio_May2014/Aminooxybiotin/P5767_P5768/201400324_Isobar_P5767.P5768.csv", header=TRUE, stringsAsFactors=F)
neg <- read.csv2("~/MP_RStudio_May2014/Aminooxybiotin/P5772_P5774/20140324_Isobar_P5772.P5774.csv", header=TRUE, stringsAsFactors=F)
pos$protein.acs <-  gsub("[,-].*",'',pos$protein.acs)   # Labelled sample contains 465 proteins
neg$protein.acs <-  gsub("[,-].*",'',neg$protein.acs)   # Ctrl contains 553 proteins
dim(pos); dim(neg)

# Show the overlap graphically
library(gplots)
pdf("P5767-8_P5772-4_Overlap.pdf")
venn(list(pos=pos$protein.acs,
          neg=neg$protein.acs))
dev.off();

# Create a combined dataframe
sc.df <- merge(pos[,c("protein.acs", "ID", "Description", "Gene", "global.spectra.count")], 
               neg[, c("ID", "global.spectra.count")], 
               by="ID")
names(sc.df) <- c('ID', "AC", "Description", "Gene", "GSC.pos", "GSC.neg")
sc.df[["AC"]] <- gsub("[,-].*",'',sc.df[['AC']])    # Combined dataframe contains 393 proteins

# Search for proteins that are detected with same GSC in positiv, labeled sample and negative ctrl
stable <- sc.df[which(sc.df$GSC.pos==sc.df$GSC.neg), ]
sum.stable.pos <- sum(stable$GSC.pos)
sum.stable.neg <- sum(stable$GSC.neg)

# Normalize the data to the sum of GSCs for stable subset
GSC.pos.norm <- sc.df$GSC.pos/sum.stable.pos
GSC.neg.norm <- sc.df$GSC.neg/sum.stable.neg
sc.df <- cbind(sc.df, GSC.pos.norm, GSC.neg.norm)

# Calculate the ratio between sample and neg. ctrl and add them to the dataframe
ratio <- sc.df$GSC.pos.norm/sc.df$GSC.neg.norm
sqrt.ratio <- sqrt(ratio)
sc.df <- cbind(sc.df, ratio, sqrt.ratio)

# Clean the data frame of NA values
sc.df <- na.omit(sc.df)
dim(sc.df) # 385 proteins left

# Look for values that are infinite
which(sc.df$ratio=="Inf")
which(sc.df$sqrt.ratio=="Inf")
# Line 156 needs to be removed
sc.df <- sc.df[-156,]
write.csv2(sc.df, file="P5767-8_P5772-4.csv")


# File created above can be read in directly ==========================
sc.df <- read.csv2("P5767-8_P5772-4.csv", header=T, stringsAsFactors=F)
dim(sc.df); str(sc.df)

# Plot density
par(mfrow=c(1,1))
plot(density(sc.df$ratio), xlab="Normalized spectral count ratios", main="Aminooxybiotin")

# Read in the PM proteins identified in the sample
PM_5767_5768 <- read.csv2("~/MP_RStudio_May2014/Aminooxybiotin/P5767_P5768/PM_P5767.P5768.csv", header=TRUE, stringsAsFactors=F)

# Identify "contaminants" in sample, i.e. proteins whose ratio is smaller than defined threshold ===============
# Take different cutoffs and plot distribution
pdf("Distrib_scRatioDiffCutoffs_sqrt_P5767-8_P5772-4.pdf")
par(mfrow=c(2,3))
thres <- c(1,1.5,2,3,4,5)
for(thre in thres){
  cont <- sc.df[which(sc.df$ratio<=thre), ] 
  clean.sf <- sc.df[-which(sc.df$ratio<=thre), ]
  extr.PM <- clean.sf[grep(paste(PM_5767_5768$protein.acs, collapse="|"), clean.sf$AC), ]
  non.extr.PM <- clean.sf[grep(paste(PM_5767_5768$protein.acs, collapse="|"), clean.sf$AC, invert=TRUE), ]
  extr.PM <- na.omit(extr.PM)
  non.extr.PM <- na.omit(non.extr.PM)
  dens.PM <- density(extr.PM$sqrt.ratio)
  plot(density(non.extr.PM$sqrt.ratio), xlim=c(0,7), ylim=c(0,3),
       xlab="Square root of normalized spectral count ratios", 
       main=paste("ratio cutoff:", thre, sep=" "))
  lines(dens.PM, col='red')
  legend("topright", legend = c('PM', 'non-PM'), fill = c("red", "black"))
}
dev.off()

# Use cutoff <=1 to define a new, clean subset
cont <- sc.df[which(sc.df$ratio<=1), ]; dim(cont)   # 207 contaminants
clean.sf <- sc.df[-which(sc.df$ratio<=1), ]
dim(sc.df); dim(cont); dim(clean.sf)

# Delete the contaminants from the sample dataset
length(intersect(cont$AC, pos$protein.acs))
pos.clean <- pos[grep(paste(cont$AC, collapse="|"), pos$protein.acs, invert=TRUE), ]; dim(pos.clean); dim(pos)
write.csv2(pos.clean, file="/slow/home/khoermann/MP/Aminooxybiotin/P5767_P5768/P5767_P5768_cleaned.csv")

# Graphical representations ==========================================================================
# For the scatterplot, it's better to leave the contaminants in
extr.PM.raw <- sc.df[grep(paste(PM_5767_5768$protein.acs, collapse="|"), sc.df$AC), ]
non.extr.PM.raw <- sc.df[grep(paste(PM_5767_5768$protein.acs, collapse="|"), sc.df$AC, invert=TRUE), ]
extr.PM.raw <- na.omit(extr.PM.raw)
non.extr.PM.raw <- na.omit(non.extr.PM.raw)
dim(extr.PM.raw); dim(non.extr.PM.raw)

par(mfrow=c(1,1))
pdf("sc_scatterplot_P5767-8_P5772-4.pdf")
plot(log10(extr.PM.raw$GSC.neg.norm), log10(extr.PM.raw$GSC.pos.norm), col="red", cex=.7,
     xlab="log10(sc(neg))", 
     ylab="log10(sc(pos))", 
     main="Aminooxybiotin")
points(log10(non.extr.PM.raw$GSC.neg.norm), log10(non.extr.PM.raw$GSC.pos.norm), cex=.7)
abline(0,1,lty=2)
legend("bottomright", legend = c('PM', 'non-PM'), fill = c("red", "black"))
legend("topleft", legend="y=x", lty=2)
dev.off()

# Lets you select the 3 non-PM proteins with highest sc(pos)/sc(neg) ratio and label them in the scatterplot
non.extr.PM.raw <- non.extr.PM.raw[order(non.extr.PM.raw[,10], decreasing=T), ]
# identify(non.extr.PM.raw$GSC.neg.norm, non.extr.PM.raw$GSC.pos.norm, n=3, labels=non.extr.PM.raw$AC)  
# More correct way is to choose a ratio cutoff
ind <- which(non.extr.PM.raw$ratio>3); 
extreme.ratio.non.PM <- non.extr.PM.raw[ind, ]; extreme.ratio.non.PM <- na.omit(extreme.ratio.non.PM)
write.xlsx2(extreme.ratio.non.PM, file="Extreme_Ratios_nonPM_Aminooxybiotin_cutoff3.xlsx")

png("sc_ratios_P5767_P5768_vs_P5772_P5774.png")
plot(density(extr.PM$ratio), col='red', xlim=c(0,35), xlab="Normalized spectral count ratios", main="Aminooxybiotin")
lines(density(non.extr.PM$ratio))
legend("topright", legend = c('PM', 'non-PM'), fill = c("red", "black"))
dev.off()
