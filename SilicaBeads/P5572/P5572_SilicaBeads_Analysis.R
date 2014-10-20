## Analysis of tryout with Silica bead protocol, 50 Mio. KBM7 cells

getwd(); setwd('~/MP_RStudio_May2014//SilicaBeads/')

# Read processed data
library(xlsx)

df.5572 <- read.xlsx2('20131203.SilicaBeads.P5572.results.xlsx', sheetIndex=1, header=TRUE, stringsAsFactors=F)
df.5572$Protein.Acs <- gsub("[,-].*",'',df.5572$Protein.Acs)  # 544 proteins identified

mem.P5572 <- df.5572[grep('membrane', df.5572$Go.Cc), ]   # 185 proteins

# Search results with a vector of GOCC annot.5572ations of interest
toMatch <- c('plasma membrane', 'cell surface', 'extracellular', 'cell membrane')
matches.P5572 <- df.5572[grep(paste(toMatch, collapse='|'), df.5572$Go.Cc), ]   # 119 matches
length(intersect(mem.P5572$Protein.Acs, matches.P5572$Protein.Acs)) # 91 non-PM membrane proteins

df.5572$matched_mask <- sapply( strsplit(df.5572$Go.Cc, ';', fixed = TRUE ),
                                  function ( CCs ) length( CCs ) > 0 && all( CCs %in% 'integral to membrane' ) )
P5572.itm <- subset(df.5572, matched_mask)   # No matches for this subset

# Find SLCs
length(grep("SLC", df.5572$Gene))   # 1 match
slc.silica <- df.5572[grep("SLC", df.5572$Gene), ]
write.csv2(slc.amino, file="SLCs_P5572.csv")    # Save SLCs to .csv file  

-----------------------------------------------------------------------------

x <- df.5572[['Protein.Acs']]

library(biomaRt)
unimart <- useMart("unimart",dataset="uniprot",
                   host="www.ebi.ac.uk",path="/uniprot/biomart/martservice")


# Take a look at the available attributes and filters
attr = listAttributes(unimart)
head(attr)
filters = listFilters(unimart)
head (filters)

go.terms=getBM(attributes=c("accession", "db2go_c__dm_description"),filters="accession",
               values=x,mart=unimart)


# Organize the output as dataframe
DF <- as.data.frame(go.terms)

# Improve the representation to show one protein per row
DF2 <- aggregate(DF[2], DF[-2], 
                 FUN = function(X) paste(unique(X), collapse=", "))   # 530 proteins annot.5572ated
colnames(DF2) <- c("Protein.Acs", "Go.Cc")

combi <- merge(df.5572[,1:12], DF2, by="Protein.Acs")     # Combine Isobar data and annot.5572ations into one dataframe
write.csv2(combi, file="P5572.GOterms.csv")

# Read in the file created before ==========================================
annot.5572 <- read.csv2("~/MP_RStudio_May2014/SilicaBeads/P5572.GOterms.csv", header=T, stringsAsFactors=F); dim(annot.5572)
annot.5572$norm.sc <- annot.5572$Global.Spectra.Count/sum(annot.5572$Global.Spectra.Count)   # add a column with the normalized sc

# Look for SLCs ------------------------------------------------------------
length(unique(grep("SLC", annot.5572$Gene)))   # 1 SLCs
length(unique(grep("ABC", annot.5572$Gene)))   # 0 ABCs
slc <- annot[grep("SLC", annot.5572$Gene),]
abc <- annot[grep("ABC", annot.5572$Gene),]
write.csv2(rbind(slc, abc), file="~/MP_RStudio_May2014/SilicaBeads/SLCsABCs_P5572.csv")

mem <- annot.5572[grep('membrane', annot.5572$Go.Cc), ]; dim(mem)  # 183 proteins

# PM only
PM.only <- c('plasma membrane', 'cell surface', 'cell membrane')
matches.pm.only <- annot.5572[grep(paste(PM.only, collapse='|'), annot.5572$Go.Cc), ]; dim(matches.pm.only)  # 91 matches
write.csv2(matches.pm.only, file="~/MP_RStudio_May2014/SilicaBeads/PMOnly_P5572.csv")

# Search results with a vector of GOCC annot.5572ations of interest
toMatch <- c('plasma membrane', 'cell surface', 'extracellular', 'cell membrane')
matches.pm <- annot.5572[grep(paste(toMatch, collapse='|'), annot.5572$Go.Cc), ]   # 258 matches
length(setdiff(mem$Protein.Acs, matches.pm$Protein.Acs)) # 65 non-PM membrane proteins
non_PM_mem <- annot.5572[grep(paste(setdiff(mem$Protein.Acs, matches.pm$Protein.Acs), collapse="|"), annot.5572$Protein.Acs), ]
excluded <- c(matches.pm$Protein.Acs, non_PM_mem$Protein.Acs)
other <- annot.5572[grep(paste(excluded, collapse="|"), annot.5572$Protein.Acs, invert=T), ]
dim(matches.pm); dim(non_PM_mem); dim(other); dim(matches.pm)[1]+dim(non_PM_mem)[1]+dim(other)[1]

write.csv2(matches.pm, file="PM_P5572.csv")
write.csv2(non_PM_mem, file="Non_PM_mem_P5572.csv")
write.csv2(other, file="Other_P5572.csv")

annot.5572$matched_mask <- sapply( strsplit(annot.5572$Go.Cc, ';', fixed = TRUE ),
                              function ( CCs ) length( CCs ) > 0 && all( CCs %in% 'integral to membrane' ) )
P5572.itm <- subset(annot.5572, matched_mask)   # No matches for this subset

