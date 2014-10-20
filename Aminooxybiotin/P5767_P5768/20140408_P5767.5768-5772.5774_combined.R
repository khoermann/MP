# Include negative control into the analysis

# Set your WD
setwd('/slow/home/khoermann/MP/Aminooxybiotin/'); getwd();

# Read in the data
pos <- read.csv2("201400324_Isobar_P5767.P5768.csv", header=TRUE, stringsAsFactors=F)
neg <- read.csv2("20140324_Isobar_P5772.P5774.csv", header=TRUE, stringsAsFactors=F)
pos$protein.acs <-  gsub("[,-].*",'',pos$protein.acs)
neg$protein.acs <-  gsub("[,-].*",'',neg$protein.acs)

# Create a combined dataframe
sc.df <- merge(pos[,c("protein.acs", "ID", "Description", "Gene", "global.spectra.count")], 
               neg[, c("ID", "global.spectra.count")], 
               by="ID")
names(sc.df) <- c('ID', "AC", "Description", "Gene", "GSC.pos", "GSC.neg")
sc.df[["AC"]] <- gsub("[,-].*",'',sc.df[['AC']])

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
intersect(cont$AC, pos$protein.acs)
pos.clean <- pos[grep(paste(cont$AC, collapse="|"), pos$protein.acs, invert=TRUE), ]

# Define the new "purified" dataset
sc.df <- sc.df[-which(sc.df$ratio<1), ]
plot(density(sc.df$ratio), lty=2)

# Read in the GO annotations for the sample
GO_pos <- read.csv2("20140324.Aminooxy.P5767.P5768-34.GOterms.csv", header=TRUE, stringsAsFactors=F)

GO_pos_clean <- GO_pos[grep(paste(pos.clean$protein.acs, collapse="|"), GO_pos$accession), ]

toMatch <- c('plasma membrane', 'cell surface', 'extracellular')
length(grep(paste(toMatch, collapse="|"), GO_pos_clean$go_name))

library(biomaRt)
unimart <- useMart("unimart",dataset="uniprot",
                   host="www.ebi.ac.uk",path="/uniprot/biomart/martservice")
GO.sc.df <- getBM(attributes=c("accession", "name", "go_name"), filters = 'accession', 
                      values = sc.df$AC, mart=unimart)
# Organize the output as dataframe
DF <- as.data.frame(GO.sc.df)

# Improve the representation to show one protein per row
DF2 <- aggregate(DF[3], DF[-3], 
                 FUN = function(X) paste(unique(X), collapse=", "))

length(grep(paste(toMatch, collapse="|"), DF2$go_name))
PM_clean.sc.df <- DF2[grep(paste(toMatch, collapse="|"), DF2$go_name), ]
PM.ratios <- sc.df[grep(paste(PM_clean.sc.df$accession, collapse="|"), sc.df$AC), "ratio"]
lines(density(PM.ratios), lty=1)
# Graphical representation
boxplot(non.extr.PM$ratio, extr.PM$ratio)

dens.PM <- density(extr.PM$ratio)
plot(dens.PM, col='red', xlim=c(0,2))
lines(density(non.extr.PM$ratio))

