## Deeper analysis of SilicaBeads data, 2nd round, 30 Mio. cells

setwd('/slow/home/khoermann/MP/SilicaBeads/'); getwd(); 
  
silica.df <- read.csv2('20140325_Isobar_P5744.csv', header=TRUE, stringsAsFactors=F)    # 117 proteins
head(silica.df)
str(silica.df)
  
silica.df$protein.acs <- gsub("[,-].*",'',silica.df$protein.acs)

library(biomaRt)
unimart <- useMart("unimart",dataset="uniprot",
                   host="www.ebi.ac.uk",path="/uniprot/biomart/martservice")

# Take a look at the available attributes and filters
attr = listAttributes(unimart)
head(attr)
filters = listFilters(unimart)
head (filters)

go.terms=getBM(attributes=c("accession", "db2go_c__dm_description"),filters="accession",
               values=silica.df$protein.acs,mart=unimart)

# Organize the output as dataframe
DF <- as.data.frame(go.terms)

# Improve the representation to show one protein per row
DF2 <- aggregate(DF[2], DF[-2], 
                 FUN = function(X) paste(unique(X), collapse=", "))   # 114 proteins annotated
colnames(DF2) <- c("protein.acs", "Go.Cc")

combi <- merge(silica.df, DF2, by="protein.acs")     # Combine Isobar data and annotations into one dataframe
write.csv2(combi, file="P5744.GOterms.csv")

annot.5744 <- read.csv2("~/MP_RStudio_May2014/SilicaBeads/P5744.GOterms.csv", header=T, stringsAsFactors=F); dim(annot.5744)
annot.5744$norm.sc <- annot.5744$global.spectra.count/sum(annot.5744$global.spectra.count)   # add a column with the normalized sc

# Look for SLCs
length(unique(grep("SLC", annot.5744$Gene)))   # 1 SLCs
slc.5744 <- annot.5744[grep("SLC", annot.5744$Gene),]
length(unique(grep("ABC", annot.5744$Gene)))   # 1 SLCs
abc.5744 <- annot.5744[grep("ABC", annot.5744$Gene),]
write.csv2(rbind(abc.5744,slc.5744), file="~/MP_RStudio_May2014/SilicaBeads/ABCsSLCs_P5744.csv")  

# PM only
PM.only <- c('plasma membrane', 'cell surface', 'cell membrane')
matches.pm.only <- annot.5744[grep(paste(PM.only, collapse='|'), annot.5744$Go.Cc), ]; dim(matches.pm.only)  # 91 matches
write.csv2(matches.pm.only, file="~/MP_RStudio_May2014/SilicaBeads/PMOnly_P5744.csv")

# Identify membrane proteins
mem <- annot.5744[grep('membrane', annot.5744$Go.Cc), ]; dim(mem)   # 30 proteins

# Search for PM proteins
toMatch <- c('plasma membrane', 'cell surface', 'extracellular', 'cell membrane')
matches.pm <- annot.5744[grep(paste(toMatch, collapse='|'), annot.5744$Go.Cc), ]    # Gives 81 matches
length(setdiff(mem$protein.acs, matches.pm$protein.acs)) # 2 non-PM membrane proteins
non_PM_mem <- annot.5744[grep(paste(setdiff(mem$protein.acs, matches.pm$protein.acs), collapse="|"), annot.5744$protein.acs), ]
excluded <- c(matches.pm$protein.acs, non_PM_mem$protein.acs)
other <- annot.5744[grep(paste(excluded, collapse="|"), annot.5744$protein.acs, invert=T), ]
dim(matches.pm); dim(non_PM_mem); dim(other); dim(matches.pm)[1]+dim(non_PM_mem)[1]+dim(other)[1]

write.csv2(matches.pm, file="PM_P5744.csv")
write.csv2(non_PM_mem, file="Non_PM_mem_P5744.csv")
write.csv2(other, file="Other_P5744.csv")


annot.5744$matched_mask <- sapply( strsplit(annot.5744$Go.Cc, ';', fixed = TRUE ),
                                   function ( CCs ) length( CCs ) > 0 && all( CCs %in% 'integral to membrane' ) )
P5744.itm <- subset(annot.5744, matched_mask)   # No matches for this subset

