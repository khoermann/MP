## Deeper analysis of 2nd round of Aminooxybiotin data

setwd('~/MP_RStudio_May2014/Aminooxybiotin/P5767_P5768/'); getwd(); 
  
aminooxy.df <- read.csv2('201400324_Isobar_P5767.P5768.csv', header=TRUE, stringsAsFactors=F)   # 465 proteins
head(aminooxy.df)
str(aminooxy.df)

 aminooxy.df$protein.acs <- gsub("[,-].*",'',aminooxy.df$protein.acs)
  
  library(biomaRt)
  unimart <- useMart("unimart",dataset="uniprot",
                     host="www.ebi.ac.uk",path="/uniprot/biomart/martservice")
  
  # Take a look at the available attributes and filters
  attr = listAttributes(unimart)
  head(attr)
  filters = listFilters(unimart)
  head (filters)
  
  go.terms=getBM(attributes=c("accession", "db2go_c__dm_description"),filters="accession",
                 values=aminooxy.df$protein.acs,mart=unimart)
  
  # Organize the output as dataframe
  DF <- as.data.frame(go.terms)
  
  # Improve the representation to show one protein per row
  DF2 <- aggregate(DF[2], DF[-2], 
                   FUN = function(X) paste(unique(X), collapse=", "))   # 457 proteins annotated
  colnames(DF2) <- c("protein.acs", "Go.Cc")
  
  combi <- merge(aminooxy.df, DF2, by="protein.acs")     # Combine Isobar data and annotations into one dataframe
  write.csv2(combi, file="P5767_5768.GOterms.csv")  

  annot.5767.5768 <- read.csv2("~/MP_RStudio_May2014/Aminooxybiotin//P5767_P5768/P5767_5768.GOterms.csv", header=T, stringsAsFactors=F); dim(annot.5767.5768)
  annot.5767.5768$norm.sc <- annot.5767.5768$global.spectra.count/sum(annot.5767.5768$global.spectra.count)   # add a column with the normalized sc

  # Look for SLCs
length(unique(grep("SLC", annot.5767.5768$Gene)))   # 5 SLCs
length(unique(grep("ABC", annot.5767.5768$Gene)))   # 1 ABCs

slc.5767.5768 <- annot.5767.5768[grep("SLC", annot.5767.5768$Gene),]
abc.5767.5768 <- annot.5767.5768[grep("ABC", annot.5767.5768$Gene),]
write.csv2(rbind(slc.5767.5768, abc.5767.5768), file="SLCsABCs_P5767.5768.csv") 
  
mem <- annot.5767.5768[grep('membrane', annot.5767.5768$Go.Cc), ]; dim(mem)   # 162 proteins
  
# PM only
PM.only <- c('plasma membrane', 'cell surface', 'cell membrane')
matches.pm.only <- annot.5767.5768[grep(paste(PM.only, collapse='|'), annot.5767.5768$Go.Cc), ]; dim(matches.pm.only)  # 50 matches
write.csv2(matches.pm.only, file="~/MP_RStudio_May2014/Aminooxybiotin/P5767_P5768/PMOnly_P5767-8.csv")


  # Search results with a vector of GOCC annotations of interest
  toMatch <- c('plasma membrane', 'cell surface', 'extracellular', 'cell membrane')
  matches.pm <- annot.5767.5768[grep(paste(toMatch, collapse='|'), annot.5767.5768$Go.Cc), ]    # Gives 280 matches
  length(setdiff(mem$protein.acs, matches.pm$protein.acs)) # 40 non-PM membrane proteins
  non_PM_mem <- annot.5767.5768[grep(paste(setdiff(mem$protein.acs, matches.pm$protein.acs), collapse="|"), annot.5767.5768$protein.acs), ]
  excluded <- c(matches.pm$protein.acs, non_PM_mem$protein.acs)
  other <- annot.5767.5768[grep(paste(excluded, collapse="|"), annot.5767.5768$protein.acs, invert=T), ]
  dim(matches.pm); dim(non_PM_mem); dim(other); dim(matches.pm)[1]+dim(non_PM_mem)[1]+dim(other)[1]
  
  write.csv2(matches.pm, file="PM_P5767.P5768.csv")
  write.csv2(non_PM_mem, file="Non_PM_mem_P5767.P5768.csv")
  write.csv2(other, file="Other_P5767.P5768.csv")
  
  annot.5767.5768$matched_mask <- sapply( strsplit(annot.5767.5768$Go.Cc, ';', fixed = TRUE ),
                                     function ( CCs ) length( CCs ) > 0 && all( CCs %in% 'integral to membrane' ) )
  P5767.5768.itm <- subset(annot.5767.5768, matched_mask)   # No matches for this subset
  
  
# Plot sc distribution ========================================
pdf("~/MP_RStudio_May2014/Aminooxybiotin/P5767_P5768/GSC_distrib.pdf")
plot(log10(matches.pm$global.spectra.count), 
     ylab="Log10(Global Spectra Count)", 
     xlab="Proteins",
     col="red")
points(log10(non_PM_mem$global.spectra.count), col="orange")
points(log10(other$global.spectra.count))
legend("topright", legend=c("PM", "non-PM Mem", "Other"), fill=c("red", "orange", "black"))
dev.off();

# "Classic" vertical boxplot
pdf("~/MP_RStudio_May2014/Aminooxybiotin/P5767_P5768/sc_distrib_vertical.pdf")
boxplot(log10(matches.pm$global.spectra.count), log10(non_PM_mem$global.spectra.count), log10(other$global.spectra.count), 
        horizontal=F,
        ylab="Log10(Global Spectra Counts)", 
        names=c("PM", "Mem non-PM", "other"), 
        col=c("red", "orange", "grey"), 
        main="Global Spectra Count distribution labelled sample")
dev.off()

# Horizontal version of the boxplot
pdf("~/MP_RStudio_May2014/Aminooxybiotin/P5767_P5768/sc_distrib_horizontal.pdf")
boxplot(log10(matches.pm$global.spectra.count), log10(non_PM_mem$global.spectra.count), log10(other$global.spectra.count), 
        horizontal=T,
        xlab="Log10(Global Spectra Counts)", 
        names=c("PM", "Mem non-PM", "other"), 
        col=c("red", "orange", "grey"), 
        main="Global Spectra Count distribution labelled sample")
dev.off()
  
