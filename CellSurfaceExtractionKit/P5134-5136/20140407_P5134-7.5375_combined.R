# Include negative control into the analysis

# Set your WD
setwd('/slow/home/khoermann/MP/CellSurfaceExtractionKit/'); getwd();

# Read in the data
library("xlsx")
pos <- read.xlsx2("20130705.mp.results.xlsx", header=TRUE, sheetIndex=1, stringsAsFactors=F)
neg <- read.xlsx2("20131202.NegCtrl.P5375.results.xlsx", header=TRUE, sheetIndex=1, stringsAsFactors=F)
pos$Protein.Acs <-  gsub("[,-].*",'',pos$Protein.Acs)
neg$Protein.Acs <-  gsub("[,-].*",'',neg$Protein.Acs)

# Create a combined dataframe
sc.df <- merge(pos[,c("Protein.Acs", "ID", "Description", "Gene", "Global.Spectra.Count")], 
               neg[, c("ID", "Global.Spectra.Count")], 
               by="ID")
names(sc.df) <- c('ID', "AC", "Description", "Gene", "GSC.pos", "GSC.neg")
sc.df[["AC"]] <- gsub("[,-].*",'',sc.df[['AC']])
sc.df[, 5:6] <- sapply(sc.df[,5:6], as.numeric)
sc.df <- sc.df[-11,]

# Normalize the data to the sum of GSC
GSC.pos.norm <- sc.df$GSC.pos/sum(sc.df$GSC.pos)
GSC.neg.norm <- sc.df$GSC.neg/sum(sc.df$GSC.neg)
sc.df <- cbind(sc.df, GSC.pos.norm, GSC.neg.norm)

# Calculate the ratio between sample and neg. ctrl
ratio <- sc.df$GSC.pos.norm/sc.df$GSC.neg.norm
sc.df <- cbind(sc.df, ratio)
# Clean the data frame of NA values
sc.df <- na.omit(sc.df)

# Plot density
plot(density(sc.df$ratio))

# Identify "contaminants" in sample, i.e. proteins whose ratio is <1
cont <- sc.df[which(sc.df$ratio<1), ]

# Delete the contaminants from the sample dataset
intersect(cont$AC, pos$Protein.Acs)
pos.clean <- pos[grep(paste(cont$AC, collapse="|"), pos$Protein.Acs, invert=TRUE), ]

toMatch <- c('plasma membrane', 'cell surface', 'extracellular')
length(grep(paste(toMatch, collapse="|"), pos.clean$Go.Cc))
length(which(pos.clean$Go.Cc=='integral to membrane'))

# Define the new "purified" dataset
sc.df <- sc.df[-which(sc.df$ratio<1), ]
plot(density(sc.df$ratio), lty=2)

# Read in the PM proteins identified in the sample
extr.PM <- sc.df[grep(paste(ACs.PM.amino.first, collapse="|"), sc.df$AC), ]
non.extr.PM <- sc.df[grep(paste(ACs.PM.amino.first, collapse="|"), sc.df$AC, invert=TRUE), ]
extr.PM <- na.omit(extr.PM)
non.extr.PM <- na.omit(non.extr.PM)

# Graphical representation
boxplot(non.extr.PM$ratio, extr.PM$ratio)
plot(non.extr.PM$ratio)
points(extr.PM$ratio, col='red')

dens.PM <- density(extr.PM$ratio)
plot(dens.PM, col='red', xlim=c(0,40))
lines(density(non.extr.PM$ratio))
