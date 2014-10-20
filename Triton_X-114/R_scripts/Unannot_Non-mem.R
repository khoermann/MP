# Identify non-membrane proteins and unannot.detated proteins in detergent and aqueous phase to maybe find potential new membrane proteins
require(xlsx)

# P5834/P5835, detergent phase --------------------------------------------------------------------
# Read in the complete dataset obtained
setwd('~/MP_RStudio_May2014/Triton_X-114/P5834_5835'); getwd(); 
det.df <- read.csv2('20140325_Isobar_P5834.P5835.csv', header=TRUE, stringsAsFactors=F)
head(det.df)
str(det.df); dim(det.df)
det.df$protein.acs <- gsub("[,-].*",'',det.df$protein.acs)

# Read in the annotated protein set
annot.det <- read.csv2("~/MP_RStudio_May2014/Triton_X-114/P5834_5835//P5834_5835.GOTerms.csv", header=T, stringsAsFactors=F); dim(annot.det)

# Extract all non-membrane-proteins 
non.mem.det <- annot.det[grep('membrane', annot.det$Go.Cc, invert=T), ]; dim(non.mem.det);    # 636 proteins
write.csv2(non.mem.det, file="~/MP_RStudio_May2014/Triton_X-114/P5834_5835/Non-mem_P5834.P5835.csv")

# Find unannotated proteins
setdiff <- setdiff(det.df$protein.acs, annot.det$protein.acs)
unannot.det <- det.df[grep(paste(setdiff, collapse="|"), det.df$protein.acs), ]
unannot.det$Go.Cc <- rep("NA", 57)

# Bring unannotated proteins and non-membrane proteins together in one file
non.mem.det <- non.mem.det[,-2]
unannot.det <- unannot.det[,-1]
unannot.det.non.mem.det <- rbind(non.mem.det, unannot.det)
write.csv2(unannot.det.non.mem.det, file="NotAnnotated_Non-mem_P5834_P5835.csv")


# P5836/P5837, aqueous phase --------------------------------------------------------------------
# Read in the complete dataset obtained
setwd('~/MP_RStudio_May2014/Triton_X-114/P5836_5837'); getwd(); 
aqu.df <- read.csv2('20140325_Isobar_P5836.P5837.csv', header=TRUE, stringsAsFactors=F)
head(aqu.df)
str(aqu.df); dim(aqu.df)
aqu.df$protein.acs <- gsub("[,-].*",'',aqu.df$protein.acs)

# Read in the annotated protein set
annot.aqu <- read.csv2("~/MP_RStudio_May2014/Triton_X-114/P5836_5837//P5836_5837_GOTerms.csv", header=T, stringsAsFactors=F); 
dim(annot.aqu)

# Extract all non-membrane-proteins 
non.mem.aqu <- annot.aqu[grep('membrane', annot.aqu$Go.Cc, invert=T), ]; dim(non.mem.aqu);    # 1201 proteins
write.csv2(non.mem.aqu, file="~/MP_RStudio_May2014/Triton_X-114/P5836_5837/Non-mem_P5836.P5837.csv")

# Find unannotated proteins
setdiff <- setdiff(aqu.df$protein.acs, annot.aqu$protein.acs)
unannot.aqu <- aqu.df[grep(paste(setdiff, collapse="|"), aqu.df$protein.acs), ]
unannot.aqu$Go.Cc <- rep("NA", 60)

# Bring unannotated proteins and non-membrane proteins together in one file
non.mem.aqu <- non.mem.aqu[,-c(1,3,4)]
unannot.aqu <- unannot.aqu[,-c(1,2)]
unannot.aqu.non.mem.aqu <- rbind(non.mem.aqu, unannot.aqu)
write.csv2(unannot.aqu.non.mem.aqu, file="NotAnnotated_Non-mem_P5836_P5837.csv")


###################################################################################################
# Continue working from here with the created files

pot.det <- read.csv2("~/MP_RStudio_May2014/Triton_X-114/P5834_5835/NotAnnotated_Non-mem_P5834_P5835.csv", header=T, stringsAsFactors=F)
pot.aqu <- read.csv2("~/MP_RStudio_May2014/Triton_X-114/P5836_5837/NotAnnotated_Non-mem_P5836_P5837.csv", header=T, stringsAsFactors=F)

# ATTENTION: There are two types of proteins that are potentially interesting: 
# Those only found in detergent phase and those found in both phases, but enriched in detergent phase

pot.det.only <- setdiff(pot.det$protein.acs, pot.aqu$protein.acs)
det.only <- pot.det[grep(paste(pot.det.only, collapse="|"), pot.det$protein.acs), ]
setwd("~/MP_RStudio_May2014/Triton_X-114/")
write.xlsx2(det.only, file="NotAnnotated_Non-mem_DetOnly_P5834.P5835.xlsx", sheetName="Proteins only found in Detergent phase")

int <- intersect(pot.det$protein.acs, pot.aqu$protein.acs)
both.det <- pot.det[grep(paste(int, collapse="|"), pot.det$protein.acs), ]
both.aqu <- pot.aqu[grep(paste(int, collapse="|"), pot.aqu$protein.acs), ]

# Create a combined dataframe
both.df <- merge(both.det[,c("protein.acs", "ID", "Description", "Gene", "global.spectra.count", "Go.Cc")], 
               both.aqu[, c("protein.acs", "global.spectra.count")], 
               by="protein.acs")
names(both.df) <- c("AC", "ID", "Description", "Gene", "GSC.det", "Go.Cc", "GSC.aqu")
both.df <- both.df[, c(1,2,3,4,6,5,7)]    # Change the order of columns to look better and have the GSC columns next to each other

# Search for proteins that are detected with same GSC in positiv, labeled sample and negative ctrl
stable <- both.df[which(both.df$GSC.det==both.df$GSC.aqu), ]
sum.stable.pos <- sum(stable$GSC.det)
sum.stable.neg <- sum(stable$GSC.aqu)

# Normalize the data to the sum of GSCs for stable subset
GSC.det.norm <- both.df$GSC.det/sum.stable.pos
GSC.aqu.norm <- both.df$GSC.aqu/sum.stable.neg
both.df <- cbind(both.df, GSC.det.norm, GSC.aqu.norm)

# Calculate the ratio between sample and neg. ctrl
ratio <- both.df$GSC.det.norm/both.df$GSC.aqu.norm
sqrt.ratio <- sqrt(ratio)
both.df <- cbind(both.df, ratio, sqrt.ratio)

# Save this to a file
write.csv2(both.df, file="~/MP_RStudio_May2014/Triton_X-114/combNotAnnot_Non-mem_P5834-P5837.csv")

# Identify those proteins that have a higher sc for detergent phase than for aqueous phase
cand <- both.df[which(both.df$ratio>1), ]
setwd("~/MP_RStudio_May2014/Triton_X-114/")
write.xlsx2(cand, file="PotMemProt_RatioCutoff1.xlsx", sheetName="PotMemProt")


# You can also start from here by reading in the data from the .csv file =================================================
both.df <- read.csv2("~/MP_RStudio_May2014/Triton_X-114//combNotAnnot_Non-mem_P5834-P5837.csv", header=T, stringsAsFactors=F)
dim(both.df)

# PLOTTING PART ####################################################################
# For plotting, the Inf ratio and NA make problems
which(both.df$ratio=="Inf")

plot.both <- both.df [-110,]
ratios.to.plot <- na.omit(plot.both$ratio)
sqrt.to.plot <- na.omit(plot.both$sqrt.ratio)
fourth.to.plot <- sqrt(sqrt.to.plot)

summary(ratios.to.plot)

# Identify "contaminants" in sample, i.e. proteins whose ratio is smaller than defined threshold
# Take different cutoffs and plot distribution
pdf(".pdf")
par(mfrow=c(2,3))
thres <- c(1,1.5,2,3,4,5)
for(thre in thres){
  cont <- both.df[which(both.df$ratio<=thre), ] 
  clean.sf <- both.df[-which(both.df$ratio<=thre), ]
  to.plot <- na.omit(clean.sf$sqrt.ratio)
  plot(density(to.plot), 
       xlab="Square root of normalized spectral count ratios", 
       main=paste("ratio cutoff:", thre, sep=" "))
    }
dev.off()

pdf("~/MP_RStudio_May2014/Triton_X-114/scDistrib.NotAnnot_Non-mem_P5834-P5837.pdf")
par(mfrow=c(1,1))
plot(log10(both.df$GSC.aqu.norm), log10(both.df$GSC.det.norm), cex=.8,
     xlab="log10(norm.sc aqueous phase)",
     ylab="log10(norm.sc detergent phase)", 
     xlim=c(-2, 1.2),
     ylim=c(-2, 1.2),
     abline(0,1, lty=2))
legend("bottomright", legend=c("y=x"), lty=2)
dev.off();


