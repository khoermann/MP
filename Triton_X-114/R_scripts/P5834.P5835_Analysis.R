## Deeper analysis of Triton X-114 data
# 30 mio. cells, detergent phase

setwd('~/MP_RStudio_May2014/Triton_X-114/P5834_5835'); getwd(); 
  
triton.df <- read.csv2('20140325_Isobar_P5834.P5835.csv', header=TRUE, stringsAsFactors=F)
head(triton.df)
str(triton.df); dim(triton.df)
triton.df$protein.acs <- gsub("[,-].*",'',triton.df$protein.acs)

# Retrieve annotation from biomaRt ==========================================
library(biomaRt)
unimart <- useMart("unimart",dataset="uniprot",
                   host="www.ebi.ac.uk",path="/uniprot/biomart/martservice")

# Take a look at the available attributes and filters
attr = listAttributes(unimart)
head(attr)
filters = listFilters(unimart)
head (filters)

go.terms=getBM(attributes=c("accession", "db2go_c__dm_description"),filters="accession",
               values=triton.df$protein.acs,mart=unimart)

# Organize the output as dataframe
DF <- as.data.frame(go.terms)

# Improve the representation to show one protein per row
DF2 <- aggregate(DF[2], DF[-2], 
                 FUN = function(X) paste(unique(X), collapse=", "))   # 1881 proteins annotated
colnames(DF2) <- c("protein.acs", "Go.Cc")

combi <- merge(triton.df, DF2, by="protein.acs")     # Combine Isobar data and annotations into one dataframe
write.csv2(combi, file="~/MP_RStudio_May2014/Triton_X-114/P5834_5835.GOTerms.csv")  # Save them to a .csv file


# Continue working with the file created above =========================================================================
annot <- read.csv2("~/MP_RStudio_May2014/Triton_X-114/P5834_5835//P5834_5835.GOTerms.csv", header=T, stringsAsFactors=F); dim(annot)
  
# Look for SLCs and ABCs
length(unique(grep("SLC", annot$Gene)))   # 71 SLCs
length(unique(grep("ABC", annot$Gene)))   # 11 ABCs
slc.5834.5835 <- annot[grep("SLC", annot$Gene),]
abc.5834.5835 <- annot[grep("ABC", annot$Gene),]
write.csv2(rbind(slc.5834.5835, abc.5834.5835), file="SLCsABCs_P5834.P5835.csv")  

# Look for subcellular annotations
mem <- annot[grep('membrane', annot$Go.Cc), ]; dim(mem);    # 1245 proteins
message("Percentage:"); round(dim(mem)[1]/dim(triton.df)[1], digits=2)    # 64%
write.csv2(mem, file='Mem_P5834.P5835.csv')

# Extract all non-membrane-proteins 
non.mem <- annot[grep('membrane', annot$Go.Cc, invert=T), ]; dim(non.mem);    # 1245 proteins
write.csv2(non.mem, file="~/MP_RStudio_May2014/Triton_X-114/P5834_5835/Non-mem_P5834.P5835.csv")

# PM only
PM.only <- c('plasma membrane', 'cell surface', 'cell membrane')
matches.pm.only <- annot[grep(paste(PM.only, collapse='|'), annot$Go.Cc), ]; dim(matches.pm.only)  # 474 matches

# Search results with a vector of GOCC annotations of interest
toMatch <- c('plasma membrane', 'cell surface', 'extracellular', 'cell membrane')
matches.pm <- annot[grep(paste(toMatch, collapse='|'), annot$Go.Cc), ]; dim(matches.pm);  # 761 proteins
message("Percentage:"); round(dim(matches.pm)[1]/dim(triton.df)[1], digits=2) # 39%

length(setdiff(mem$protein.acs, matches.pm$protein.acs)) # 670 non-PM membrane proteins
non_PM_mem <- annot[grep(paste(setdiff(mem$protein.acs, matches.pm$protein.acs), collapse="|"), annot$protein.acs), ]
excluded <- c(matches.pm$protein.acs, non_PM_mem$protein.acs)
other <- annot[grep(paste(excluded, collapse="|"), annot$protein.acs, invert=T), ]
dim(matches.pm); dim(non_PM_mem); dim(other); dim(matches.pm)[1]+dim(non_PM_mem)[1]+dim(other)[1]

write.csv2(matches.pm.only, file="PMOnly_P5834.P5835.csv")
write.csv2(matches.pm, file='PM_P5834.P5835.csv')  
write.csv2(non_PM_mem, file="Non_PM_mem_P5834.P5835.csv")
write.csv2(other, file="Other_P5834.P5835.csv")

annot$matched_mask <- sapply( strsplit(annot$Go.Cc, ';', fixed = TRUE ),
                              function ( CCs ) length( CCs ) > 0 && all( CCs %in% 'integral to membrane' ) )
triton.itm <- subset(annot, matched_mask); dim(triton.itm)  # No matches for this subset

