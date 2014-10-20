## Deeper analysis of Aminooxybiotin data

  setwd('/slow/home/khoermann/MP/Aminooxybiotin/P5396_P5455/'); getwd(); 
  library(xlsx)

aminooxy.df <- read.xlsx2('20131202.Aminooxy.P5396.P5455.results.xlsx', sheetIndex=1, stringsAsFactors=F)   # 318 proteins
aminooxy.df$Protein.Acs <- gsub("[,-].*",'',aminooxy.df$Protein.Acs)
  
# Search results with a vector of GOCC annotations of interest
toMatch <- c('plasma membrane', 'cell surface', 'extracellular', 'cell membrane')
matches.pm <- aminooxy.df[grep(paste(toMatch, collapse='|'), aminooxy.df$Go.Cc), ]  # 193 matches

  # Extract those proteins with only one Go.CC: "integral to membrane"
  aminooxy.df$matched_mask <- sapply( strsplit( as.character( aminooxy.df$Go.Cc ), ';', fixed = TRUE ),
                                    function ( CCs ) length( CCs ) > 0 && all( CCs %in% 'integral to membrane' ) )
  aminooxy.itm <- subset(aminooxy.df, matched_mask)   # 7 additional ITMs

PM_aminooxy <- rbind(matches.pm, aminooxy.itm)    # Create a complete dataset (200 proteins)
write.csv2(PM_aminooxy, file="PM_P5396.P5455.csv")  # Save it to a .csv file
  
matches.mem <- aminooxy.df[grep("membrane", aminooxy.df$Go.Cc), ]   # 204 membrane
length(setdiff(matches.mem$Protein.Acs, PM_aminooxy$Protein.Acs))   # Only 21 non-PM membrane proteins

  # Identify non-PM proteins  
matches.not.pm <- aminooxy.df[grep(paste(toMatch, collapse='|'), aminooxy.df$Go.Cc, invert=TRUE), ]
  


  #####################################
  # Get annotations from biomaRt
  #####################################
  
  x <- aminooxy.df[['Protein.Acs']]
  
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
                   FUN = function(X) paste(unique(X), collapse=", "))   # 316 proteins annotated
  colnames(DF2) <- c("Protein.Acs", "Go.Cc")
  
  annot <- merge(aminooxy.df[,1:15], DF2, by="Protein.Acs")     # annotne Isobar data and annotations into one dataframe
  write.csv2(annot, file="P5396_P5455.GOterms.csv")
  
  annot <- read.csv2("~/MP_RStudio_May2014/Aminooxybiotin/P5396_P5455/P5396_P5455.GOterms.csv", header=T, stringsAsFactors=F)
  annot$norm.sc <- annot$Global.Spectra.Count/sum(annot$Global.Spectra.Count)   # add a column with the normalized sc
  
  # Look for SLCs ------------------------------------------------------------
  length(unique(grep("SLC", annot$Gene)))   # 29 SLCs
  length(unique(grep("ABC", annot$Gene)))   # 2 ABCs
  slc <- annot[grep("SLC", annot$Gene),]
  abc <- annot[grep("ABC", annot$Gene),]
  write.csv2(rbind(slc,abc), file="~/MP_RStudio_May2014/Aminooxybiotin//P5396_P5455/ABCsSLCs_P5396.P5455.csv")    # Save to .csv file
  
  # Membrane
  mem <- annot[grep('membrane', annot$Go.Cc), ]; dim(mem)  # 205 proteins
  
    # PM only
  PM.only <- c('plasma membrane', 'cell surface', 'cell membrane')
  matches.pm.only <- annot[grep(paste(PM.only, collapse='|'), annot$Go.Cc), ]; dim(matches.pm.only)  # 173 matches
  write.csv2(matches.pm.only, file="~/MP_RStudio_May2014/Aminooxybiotin//P5396_P5455/PMOnly_P5396_5455.csv")
  
  # Search results with a vector of GOCC annotations of interest
  toMatch <- c('plasma membrane', 'cell surface', 'extracellular', 'cell membrane')
  matches.pm <- annot[grep(paste(toMatch, collapse='|'), annot$Go.Cc), ]; dim(matches.pm)   # 256 matches
  length(setdiff(mem$Protein.Acs, matches.pm$Protein.Acs)) # 17 non-PM membrane proteins
  non_PM_mem <- annot[grep(paste(setdiff(mem$Protein.Acs, matches.pm$Protein.Acs), collapse="|"), annot$Protein.Acs), ]
  excluded <- c(matches.pm$Protein.Acs, non_PM_mem$Protein.Acs)
  other <- annot[grep(paste(excluded, collapse="|"), annot$Protein.Acs, invert=T), ]
  dim(matches.pm); dim(non_PM_mem); dim(other); dim(matches.pm)[1]+dim(non_PM_mem)[1]+dim(other)[1]
  
  annot$matched_mask <- sapply( strsplit(annot$Go.Cc, ';', fixed = TRUE ),
                                function ( CCs ) length( CCs ) > 0 && all( CCs %in% 'integral to membrane' ) )
  aminooxy.itm <- subset(annot, matched_mask)   # No matches for this subset
  
  write.csv2(matches.pm, file="PM_5396.P5455.csv")
  write.csv2(non_PM_mem, file="Non_PM_mem_P5396.P5455.csv")
  write.csv2(other, file="Other_P5396.P5455.csv")
-----------------------------------------------------------------------  
  # Identify 'contaminations' by other membrane-bound organelles
  length(grep(paste(c('nuclear', 'nucleus'), collapse='|'), matches.not.pm$Go.Cc))   # absolute
  length(grep(paste(c('nuclear', 'nucleus'), collapse='|'), matches.not.pm$Go.Cc))/dim(aminooxy.df)[1]    # relative
  post.nuclear <- matches.not.pm[grep(paste(c('nuclear', 'nucleus'), collapse='|'), matches.not.pm$Go.Cc, invert=TRUE), ]
  
  length(grep(paste(c('cytosol', 'cytoplasm'), collapse='|'), post.nuclear$Go.Cc))   # absolute
  length(grep(paste(c('cytosol', 'cytoplasm'), collapse='|'), post.nuclear$Go.Cc))/dim(aminooxy.df)[1]    # relative
  post.nuclear.cyto <- post.nuclear[grep(paste(c('cytosol', 'cytoplasm'), collapse='|'), post.nuclear$Go.Cc, invert=TRUE), ]
  
  length(grep("endoplasmic reticulum", post.nuclear.cyto$Go.Cc))   # absolute
  length(grep("endoplasmic reticulum", post.nuclear.cyto$Go.Cc))/dim(aminooxy.df)[1]    # relative
  
  length(grep("Golgi", post.nuclear.cyto$Go.Cc))   # absolute
  length(grep("Golgi", post.nuclear.cyto$Go.Cc))/dim(aminooxy.df)[1]    # relative
  
  length(grep(paste(c("mitochondria","mitochondrial"), collapse="|"), post.nuclear.cyto$Go.Cc))   # absolute
  length(grep(paste(c("mitochondria","mitochondrial"), collapse="|"), post.nuclear.cyto$Go.Cc))/dim(aminooxy.df)[1]   # absolute
  

# Extract Global Spectra Counts for the two groups 
  # Turn Global Spectra Counts into a numeric vector
  gsc.pm <- as.numeric(as.character(matches.pm$Global.Spectra.Count))
  gsc.not.pm <- as.numeric(as.character(matches.not.pm$Global.Spectra.Count))
  
  # Plot the data
  par(mfrow=c(1,1))
  plot(gsc.not.pm, cex=0.75, ylab='Global Spectra Count')
  points(gsc.pm, col='red', cex=0.75)
  legend(locator(1), c('Non-PM', "PM"), fill=c('black', 'red'), cex=rep(0.75,2))
  
# Display data distribution graphically
  par(mfrow=c(1,1))
  boxplot(gsc.pm, gsc.not.pm, ylab='Global Spectra Count', names=c('PM', 'Non-PM'))

par(mfrow=c(1,2))  
hist(gsc.pm,
     xlim=c(0, 500),
     xlab='Global Spectra Count',
     main='Histogram of PM protein distribution')
hist(gsc.not.pm,
     ylim=c(0, 150),
     xlab='Global Spectra Count',
     main='Histogram of non-PM protein distribution')
   
# Get parameters for the two groups
summary(gsc.pm)
summary(gsc.not.pm)
