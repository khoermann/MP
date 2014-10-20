## Deeper analysis of Pierce kit data
# Modification: pre deglycosylation by PNGaseF

setwd('/slow/home/khoermann/MP/CellSurfaceExtractionKit/P5824/'); getwd(); 
  

kit.df <- read.csv2('20140325_Isobar_P5824.csv', header=TRUE, stringsAsFactors=F)
head(kit.df)
str(kit.df)
kit.df$protein.acs <- gsub("[,-].*",'',kit.df$protein.acs)

library(biomaRt)
unimart <- useMart("unimart",dataset="uniprot",
                   host="www.ebi.ac.uk",path="/uniprot/biomart/martservice")

# Take a look at the available attributes and filters
attr = listAttributes(unimart)
head(attr)
filters = listFilters(unimart)
head (filters)

go.terms=getBM(attributes=c("accession", "db2go_c__dm_description"),filters="accession",
               values=kit.df$protein.acs,mart=unimart)

# Organize the output as dataframe
DF <- as.data.frame(go.terms)

# Improve the representation to show one protein per row
DF2 <- aggregate(DF[2], DF[-2], 
                 FUN = function(X) paste(unique(X), collapse=", "))   # 783 proteins annotated
colnames(DF2) <- c("protein.acs", "Go.Cc")

combi <- merge(kit.df, DF2, by="protein.acs")     # Combine Isobar data and annotations into one dataframe
write.csv2(combi, file="P5824.GOterms.csv")

# Start from here with GO Annotation file ======================================
annot.5824 <- read.csv2("~/MP_RStudio_May2014/CellSurfaceExtractionKit//P5824/P5824.GOterms.csv", header=T, stringsAsFactors=F)
annot.5824$norm.sc <- annot.5824$global.spectra.count/sum(annot.5824$global.spectra.count)   # add a column with the normalized sc

# Look for SLCs
length(unique(grep("SLC", annot.5824$Gene)))   # 10 SLCs
length(unique(grep("ABC", annot.5824$Gene)))   # 0 ABCs

slc.5824 <- annot.5824[grep("SLC", annot.5824$Gene),]
# abc.5824 <- annot.5824[grep("ABC", annot.5824$Gene),]
write.csv2(slc.5824, file="SLCs_P5824.csv")  

mem <- annot.5824[grep('membrane', annot.5824$Go.Cc), ]; dim(mem)  # 380 proteins

# PM only
PM.only <- c('plasma membrane', 'cell surface', 'cell membrane')
matches.pm.only <- annot.5824[grep(paste(PM.only, collapse='|'), annot.5824$Go.Cc), ]; dim(matches.pm.only)  # 165 matches

# Search results with a vector of GOCC annotations of interest
toMatch <- c('plasma membrane', 'cell surface', 'extracellular', 'cell membrane')
matches.pm <- annot.5824[grep(paste(toMatch, collapse='|'), annot.5824$Go.Cc), ]; dim(matches.pm)  # Gives 391 matches
length(setdiff(mem$protein.acs, matches.pm$protein.acs)) # 163 non-PM membrane proteins
non_PM_mem <- annot.5824[grep(paste(setdiff(mem$protein.acs, matches.pm$protein.acs), collapse="|"), annot.5824$protein.acs), ]
excluded <- c(matches.pm$protein.acs, non_PM_mem$protein.acs)
other <- annot.5824[grep(paste(excluded, collapse="|"), annot.5824$protein.acs, invert=T), ]
dim(matches.pm); dim(non_PM_mem); dim(other); dim(matches.pm)[1]+dim(non_PM_mem)[1]+dim(other)[1]

annot.5824$matched_mask <- sapply( strsplit(annot.5824$Go.Cc, ';', fixed = TRUE ),
                                   function ( CCs ) length( CCs ) > 0 && all( CCs %in% 'integral to membrane' ) )
P5824.itm <- subset(annot.5824, matched_mask)   # No matches for this subset

setwd('~/MP_RStudio_May2014//CellSurfaceExtractionKit/P5824/'); getwd(); 
write.csv(matches.pm.only, file="PMOnly_P5824.csv")
write.csv2(matches.pm, file="PM_5824.csv")
write.csv2(non_PM_mem, file="Non_PM_mem_P5824.csv")
write.csv2(other, file="Other_P5824.csv")  
  
