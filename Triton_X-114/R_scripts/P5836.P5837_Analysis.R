## Deeper analysis of Triton X-114 data
# 30 Mio. cells, acqueous phase, P5836-7

setwd('~/MP_RStudio_May2014/Triton_X-114/P5836_5837'); getwd(); 
  
acq.df <- read.csv2('20140325_Isobar_P5836.P5837.csv', header=TRUE, stringsAsFactors=F)
head(acq.df)
str(acq.df); dim(acq.df)
  
acq.df$protein.acs <- gsub("[,-].*",'',acq.df$protein.acs)

library(biomaRt)
unimart <- useMart("unimart",dataset="uniprot",
                   host="www.ebi.ac.uk",path="/uniprot/biomart/martservice")

# Take a look at the available attributes and filters
attr = listAttributes(unimart)
head(attr)
filters = listFilters(unimart)
head (filters)

go.terms=getBM(attributes=c("accession", "db2go_c__dm_description"),filters="accession",
               values=acq.df$protein.acs,mart=unimart)

# Organize the output as dataframe
DF <- as.data.frame(go.terms)

# Improve the representation to show one protein per row
DF2 <- aggregate(DF[2], DF[-2], 
                 FUN = function(X) paste(unique(X), collapse=", "))   # 1775 proteins annotated
colnames(DF2) <- c("protein.acs", "Go.Cc"); dim(DF2)

combi <- merge(acq.df, DF2, by="protein.acs"); dim(combi)     # Combine Isobar data and annotations into one dataframe
write.csv2(combi, file="P5836_5837_GOTerms.csv")


# Continue working with the file created above =========================================
annot <- read.csv2("P5836_5837_GOTerms.csv", header=T, stringsAsFactors=F); dim(annot)

# Look for SLCs and ABCs
length(unique(grep("SLC", annot$Gene)))   # 9 SLCs
length(unique(grep("ABC", annot$Gene)))   # 2 ABCs
slc.5836.5837 <- annot[grep("SLC", annot$Gene),]
abc.5836.5837 <- annot[grep("ABC", annot$Gene),]
write.csv2(rbind(slc.5836.5837, abc.5836.5837), file="SLCsABCs_P5836.P5837.csv")  

# Create "subdataframes" and files containing different subcellular annotation ==========
# Total Membrane fraction
mem <- annot[grep('membrane', annot$Go.Cc), ]; dim(mem);    # 585 proteins
message("Percentage:"); round(dim(mem)[1]/dim(acq.df)[1], digits=2)    # 32%
write.csv2(mem, file='Mem_P5836.P5837.csv')

# PM only
PM.only <- c('plasma membrane', 'cell surface', 'cell membrane')
matches.pm.only <- annot[grep(paste(PM.only, collapse='|'), annot$Go.Cc), ]; dim(matches.pm.only)  # 269 matches

# Search results with a vector of GOCC annotations of interest
toMatch <- c('plasma membrane', 'cell surface', 'extracellular', 'cell membrane')
matches.pm <- annot[grep(paste(toMatch, collapse='|'), annot$Go.Cc), ]; dim(matches.pm);  # 696 proteins
message("Percentage:"); round(dim(matches.pm)[1]/dim(acq.df)[1], digits=2) # 38%

length(setdiff(mem$protein.acs, matches.pm$protein.acs)) # 237 non-PM membrane proteins
non_PM_mem <- annot[grep(paste(setdiff(mem$protein.acs, matches.pm$protein.acs), collapse="|"), annot$protein.acs), ]
excluded <- c(matches.pm$protein.acs, non_PM_mem$protein.acs)
other <- annot[grep(paste(excluded, collapse="|"), annot$protein.acs, invert=T), ]
dim(matches.pm); dim(non_PM_mem); dim(other); dim(matches.pm)[1]+dim(non_PM_mem)[1]+dim(other)[1]

write.csv2(matches.pm.only, file="PMOnly_P5836.P5837.csv")
write.csv2(matches.pm, file='PM_P5836.P5837.csv')  
write.csv2(non_PM_mem, file="Non_PM_mem_P5836.P5837.csv")
write.csv2(other, file="Other_P5836.P5837.csv")

annot$matched_mask <- sapply( strsplit(annot$Go.Cc, ';', fixed = TRUE ),
                              function ( CCs ) length( CCs ) > 0 && all( CCs %in% 'integral to membrane' ) )
triton.itm <- subset(annot, matched_mask); dim(triton.itm)  # No matches for this subset



