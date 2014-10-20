## Deeper analysis of Triton X-114 data
# 30 Mio. cells, ctrl -> (NO phase separation)

setwd('~/MP_RStudio_May2014/Triton_X-114/P5838_5839/'); getwd(); 
  
ctrl.df <- read.csv2('20140325_Isobar_P5838.P5839.csv', header=TRUE, stringsAsFactors=F)
head(ctrl.df)
str(ctrl.df); dim(ctrl.df)  # 2036 proteins in total
  
ctrl.df$protein.acs <- gsub("[,-].*",'',ctrl.df$protein.acs)

library(biomaRt)
unimart <- useMart("unimart",dataset="uniprot",
                   host="www.ebi.ac.uk",path="/uniprot/biomart/martservice")

# Take a look at the available attributes and filters
attr = listAttributes(unimart)
head(attr)
filters = listFilters(unimart)
head (filters)

go.terms=getBM(attributes=c("accession", "db2go_c__dm_description"), filters="accession",
               values=ctrl.df$protein.acs, mart=unimart)

# Organize the output as dataframe
DF <- as.data.frame(go.terms)

# Improve the representation to show one protein per row
DF2 <- aggregate(DF[2], DF[-2], 
                 FUN = function(X) paste(unique(X), collapse=", "))   # 1969 proteins annotated
colnames(DF2) <- c("protein.acs", "Go.Cc"); dim(DF2)

combi <- merge(ctrl.df, DF2, by="protein.acs"); dim(combi)     # Combine Isobar data and annotations into one dataframe
setwd("/slow/home/khoermann/MP/Triton_X-114/");
write.csv2(combi, file="P5838_P5839.GOterms.csv")

# Use the file created before ========================================================================
annot <- read.csv2("~/MP_RStudio_May2014/Triton_X-114//P5838_5839/P5838_P5839.GOterms.csv", header=T, stringsAsFactors=F); dim(annot)
annot$norm.sc <- annot$global.spectra.count/sum(annot$global.spectra.count)   # add a column with the normalized sc

# Membrane
mem <- annot[grep('membrane', annot$Go.Cc), ]; dim(mem);    # 801 proteins
message("Percentage:"); round(dim(mem)[1]/dim(ctrl.df)[1], digits=2)    # 28%
write.csv2(mem, file='Mem_P5838.P5839.csv')

# PM only
PM.only <- c('plasma membrane', 'cell surface', 'cell membrane')
matches.pm.only <- annot[grep(paste(PM.only, collapse='|'), annot$Go.Cc), ]; dim(matches.pm.only)  # 349 matches
write.csv2(matches.pm.only, file="~/MP_RStudio_May2014/Triton_X-114//P5838_5839/PMOnly_P5838.P5839.csv")

# Search results with a vector of GOCC annotations of interest
toMatch <- c('plasma membrane', 'cell surface', 'extracellular', 'cell membrane')
matches.pm <- annot[grep(paste(toMatch, collapse='|'), annot$Go.Cc), ]; dim(matches.pm);  # 696 proteins
message("Percentage:"); round(dim(matches.pm)[1]/dim(ctrl.df)[1], digits=2) # 38%
write.csv2(matches.pm, file='PM_P5838.P5839.csv')  

length(setdiff(mem$protein.acs, matches.pm$protein.acs)) # 352 non-PM membrane proteins
non_PM_mem <- annot[grep(paste(setdiff(mem$protein.acs, matches.pm$protein.acs), collapse="|"), annot$protein.acs), ]
excluded <- c(matches.pm$protein.acs, non_PM_mem$protein.acs)
other <- annot[grep(paste(excluded, collapse="|"), annot$protein.acs, invert=T), ]
dim(matches.pm); dim(non_PM_mem); dim(other); dim(matches.pm)[1]+dim(non_PM_mem)[1]+dim(other)[1]

write.csv2(non_PM_mem, file="Non_PM_mem_P5838_5839.csv")
write.csv2(other, file="Other_P5838_5839.csv")

annot$matched_mask <- sapply( strsplit(annot$Go.Cc, ';', fixed = TRUE ),
                              function ( CCs ) length( CCs ) > 0 && all( CCs %in% 'integral to membrane' ) )
triton.itm <- subset(annot, matched_mask); dim(triton.itm)  # No matches for this subset

# Look for SLCs and ABCs
length(unique(grep("SLC", annot$Gene)))   # 16 SLCs
length(unique(grep("ABC", annot$Gene)))   # 4 ABCs
slc.5838.5839 <- annot[grep("SLC", annot$Gene),]
abc.5838.5839 <- annot[grep("ABC", annot$Gene),]
write.csv2(rbind(slc.5838.5839, abc.5838.5839), file="SLCsABCs_P5838.P5839.csv")  

