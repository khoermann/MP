## Deeper analysis of Pierce kit data
# Modification: +TX-114 phase separation, aqueous phase

setwd('~/MP_RStudio_May2014//CellSurfaceExtractionKit/P5830/'); getwd(); 

kit.df <- read.csv2('20140325_Isobar_P5830.csv', header=TRUE, stringsAsFactors=F)
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
                 FUN = function(X) paste(unique(X), collapse=", "))   # 636 proteins annotated
colnames(DF2) <- c("protein.acs", "Go.Cc")
dim(DF2)

combi <- merge(kit.df, DF2, by="protein.acs")     # Combine Isobar data and annotations into one dataframe
write.csv2(combi, file="P5830.GOterms.csv")

annot.5830 <- read.csv2("P5830.GOterms.csv", header=T, stringsAsFactors=F)
annot.5830$norm.sc <- annot.5830$global.spectra.count/sum(annot.5830$global.spectra.count)   # add a column with the normalized sc

# Look for SLCs
length(unique(grep("SLC", annot.5830$Gene)))   # 8 SLCs
length(unique(grep("ABC", annot.5830$Gene)))   # 1 ABCs

slc.5830 <- annot.5830[grep("SLC", annot.5830$Gene),]
abc.5830 <- annot.5830[grep("ABC", annot.5830$Gene),]
write.csv2(rbind(slc.5830, abc.5830), file="SLCsABCs_P5830.csv")  

mem <- annot.5830[grep('membrane', annot.5830$Go.Cc), ]; dim(mem)  # 261 proteins

# PM only
PM.only <- c('plasma membrane', 'cell surface', 'cell membrane')
matches.pm.only <- annot.5830[grep(paste(PM.only, collapse='|'), annot.5830$Go.Cc), ]; dim(matches.pm.only)  # 158 matches
write.csv2(matches.pm.only, file="~/MP_RStudio_May2014/CellSurfaceExtractionKit//P5830/PMOnly_P5830.csv")

# Search results with a vector of GOCC annotations of interest
toMatch <- c('plasma membrane', 'cell surface', 'extracellular', 'cell membrane')
matches.pm <- annot.5830[grep(paste(toMatch, collapse='|'), annot.5830$Go.Cc), ]; dim(matches.pm)  # Gives 391 matches
length(setdiff(mem$protein.acs, matches.pm$protein.acs)) # 163 non-PM membrane proteins
non_PM_mem <- annot.5830[grep(paste(setdiff(mem$protein.acs, matches.pm$protein.acs), collapse="|"), annot.5830$protein.acs), ]
excluded <- c(matches.pm$protein.acs, non_PM_mem$protein.acs)
other <- annot.5830[grep(paste(excluded, collapse="|"), annot.5830$protein.acs, invert=T), ]
dim(matches.pm); dim(non_PM_mem); dim(other); dim(matches.pm)[1]+dim(non_PM_mem)[1]+dim(other)[1]

annot.5830$matched_mask <- sapply( strsplit(annot.5830$Go.Cc, ';', fixed = TRUE ),
                                   function ( CCs ) length( CCs ) > 0 && all( CCs %in% 'integral to membrane' ) )
P5830.itm <- subset(annot.5830, matched_mask)   # No matches for this subset

write.csv2(matches.pm, file="PM_5830.csv")
write.csv2(non_PM_mem, file="Non_PM_mem_P5830.csv")
write.csv2(other, file="Other_P5830.csv")  