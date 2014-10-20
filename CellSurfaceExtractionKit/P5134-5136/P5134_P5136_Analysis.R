## Deeper analysis of data obtained with Pierce cell surface extraction kit

  setwd('/slow/home/khoermann/MP/CellSurfaceExtractionKit/P5134-5136/'); getwd(); 
  library(xlsx)

biotin.df <- read.xlsx2('20130705.mp.results.xlsx', sheetIndex=1, stringsAsFactors=F)   # Complete dataset: 2283 proteins
  biotin.df$Protein.Acs <- gsub("[,-].*",'',biotin.df$Protein.Acs)
  
# Search results with a vector of GOCC annotations of interest
toMatch <- c('plasma membrane', 'cell surface', 'extracellular', 'cell membrane')
biotin.matches.pm <- biotin.df[grep(paste(toMatch, collapse='|'), biotin.df$Go.Cc), ]   # 649 matches
biotin.matches.pm$Protein.Acs <-  gsub("[,-].*",'',biotin.matches.pm$Protein.Acs)   # Processing...
  length(unique(biotin.matches.pm$Protein.Acs)); # 561 unique ACs
    
# In addition, add those that only have one GO Annotation 'integral to membrane'
itm.kit.first <- read.xlsx2('20131203_ITM_CellSurfaceKit.xlsx', sheetIndex=1, header=TRUE, stringsAsFactors=F)  # 46 entries
acc.itm.first <- biotin.df[grep(paste(itm.kit.first$ID, collapse="|"), biotin.df$ID),2]
PM.Kit.first <- unique(c(acc.itm.first, biotin.matches.pm$Protein.Acs))   # 603 PM proteins
write.csv2(PM.Kit.first, file="PM_ACs_P5134-7.csv")  
  
biotin.mem <- biotin.df[grep("membrane", biotin.df$Go.Cc), ]    # 929 matches
biotin.mem$Protein.Acs <-  gsub("[,-].*",'',biotin.mem$Protein.Acs)   # Processsing...
  length(unique(biotin.mem$Protein.Acs))    # 807 unique ACs
length(intersect(biotin.mem$Protein.Acs, c(biotin.matches.pm$Protein.Acs, acc.itm.first)))
# Identify non-PM proteins  
biotin.matches.not.pm <- biotin.df[grep(paste(toMatch, collapse='|'), biotin.df$Go.Cc, invert=TRUE), ]




  #####################################
  # Get annotations from biomaRt
  #####################################
  
  x <- biotin.df[['Protein.Acs']]
  
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
                   FUN = function(X) paste(unique(X), collapse=", "))   # 1891 proteins annotated
  colnames(DF2) <- c("Protein.Acs", "Go.Cc")
  
  combi <- merge(biotin.df[,1:15], DF2, by="Protein.Acs")     # Combine Isobar data and annotations into one dataframe
  write.csv2(combi, file="P5134_P5136.GOterms.csv")
  
  
  # Start working with file created before ===================================================================
  annot <- read.csv2("~/MP_RStudio_May2014/CellSurfaceExtractionKit/P5134-5136/P5134_P5136.GOterms.csv", header=T, stringsAsFactors=F); dim(annot)
  annot$norm.sc <- annot$Global.Spectra.Count/sum(annot$Global.Spectra.Count)   # add a column with the normalized sc
  
  # Look for SLCs ------------------------------------------------------------
  length(unique(grep("SLC", annot$Gene)))   # 30 SLCs
  length(unique(grep("ABC", annot$Gene)))   # 14 ABCs
  slc <- annot[grep("SLC", annot$Gene),]
  abc <- annot[grep("ABC", annot$Gene),]
  write.csv2(rbind(slc, abc), file="~/MP_RStudio_May2014/CellSurfaceExtractionKit//P5134-5136/SLCsABCs_P5134-6.csv")
  
  mem <- annot[grep('membrane', annot$Go.Cc), ]   # 944 proteins
  
  # PM only
  PM.only <- c('plasma membrane', 'cell surface', 'cell membrane')
  matches.pm.only <- annot[grep(paste(PM.only, collapse='|'), annot$Go.Cc), ]; dim(matches.pm.only)  # 565 matches
  write.csv2(matches.pm.only, file="~/MP_RStudio_May2014/CellSurfaceExtractionKit//P5134-5136/PMOnly_P5134-6.csv")
  
  # Search results with a vector of GOCC annotations of interest
  toMatch <- c('plasma membrane', 'cell surface', 'extracellular', 'cell membrane')
  matches.pm <- annot[grep(paste(toMatch, collapse='|'), annot$Go.Cc), ]   # 968 matches
  length(setdiff(mem$Protein.Acs, matches.pm$Protein.Acs)) # 263 non-PM membrane proteins
  non_PM_mem <- annot[grep(paste(setdiff(mem$Protein.Acs, matches.pm$Protein.Acs), collapse="|"), annot$Protein.Acs), ]
  excluded <- c(matches.pm$Protein.Acs, non_PM_mem$Protein.Acs)
  other <- annot[grep(paste(excluded, collapse="|"), annot$Protein.Acs, invert=T), ]
  dim(matches.pm); dim(non_PM_mem); dim(other); dim(matches.pm)[1]+dim(non_PM_mem)[1]+dim(other)[1]
    
  annot$matched_mask <- sapply( strsplit(annot$Go.Cc, ';', fixed = TRUE ),
                                function ( CCs ) length( CCs ) > 0 && all( CCs %in% 'integral to membrane' ) )
  biotin.itm <- subset(annot, matched_mask)   # No matches for this subset
  
  write.csv2(matches.pm, file="PM_P5134_5136.csv")
  write.csv2(non_PM_mem, file="Non_PM_mem_P5134_5136.csv")
  write.csv2(other, file="Other_P5134_5136.csv")
----------------------------------------------------------------------------------------------
  
# Extract Global Spectra Counts for the two groups 
  # Turn Global Spectra Counts into a numeric vector
  gsc.biotin.pm <- as.numeric(as.character(biotin.matches.pm$Global.Spectra.Count))
  gsc.biotin.not.pm <- as.numeric(as.character(biotin.matches.not.pm$Global.Spectra.Count))
  
  # Plot the data
  par(mfrow=c(1,1))
  plot(gsc.biotin.not.pm, cex=0.75, ylab='Global Spectra Count')
  points(gsc.biotin.pm, col='red', cex=0.75)
  legend(locator(1), c('Non-PM', "PM"), fill=c('black', 'red'), cex=rep(0.75,2))
  
# Display data distribution graphically
  par(mfrow=c(1,1))
  boxplot(gsc.biotin.pm, gsc.biotin.not.pm, ylab='Global Spectra Count', names=c('PM', 'Non-PM'))

par(mfrow=c(1,2))  
hist(gsc.biotin.pm,
     main='Histogram of PM protein distribution',
     ylim=c(0,1500),
     xlab='Global Spectra Count')
hist(gsc.biotin.not.pm,
     main='Histogram of non-PM protein distribution', 
     xlim=c(0,400),
     xlab='Global Spectra Count')
   
# Get parameters for the two groups
summary(gsc.biotin.pm)
summary(gsc.biotin.not.pm)
  
------------------------------------
# Read in KBM7 proteome
kbm7.data <- read.csv2('/slow/home/khoermann/MP/P4569_kbm7.csv')

# Identify those proteins IDs present in both datasets
biotin.pm.ids <- as.character(biotin.matches.pm$ID)
  # Extract the PM proteins from KBM7
pm.kbm7 <- kbm7.data[grep(paste(biotin.pm.ids, collapse='|'), kbm7.data$ID), ]

# Extract sc for KBM7 
sc.pm.kbm7 <- as.numeric(pm.kbm7$Sc_mean_grp0)
  
par(mfrow=c(1,1))
boxplot(sc.pm.kbm7, gsc.biotin.pm, 
        ylab='Spectral Counts', 
        names=c('KBM7 proteome', 'Biotinylation extraction protocol'),
        main='Spectral Counts distribution for PM proteins')

summary(sc.pm.kbm7)
summary(gsc.biotin.pm)
