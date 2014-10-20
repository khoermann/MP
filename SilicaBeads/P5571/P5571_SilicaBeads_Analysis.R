## Analysis of tryout with Silica bead protocol, 30 Mio. cells

setwd('/slow/home/khoermann/MP/SilicaBeads/'); getwd();

# Read processed data
library(xlsx)

silica.df <- read.xlsx2('20131203.SilicaBeads.P5571.results.xlsx', sheetIndex=1, header=TRUE, stringsAsFactors=F)
silica.df$Protein.Acs <- gsub("[,-].*",'',silica.df$Protein.Acs); dim(silica.df)  # 277 proteins identified

mem.P5571 <- silica.df[grep('membrane', silica.df$Go.Cc), ]   # 94 proteins

# Search results with a vector of GOCC annotations of interest
toMatch <- c('plasma membrane', 'cell surface', 'extracellular', 'cell membrane')
matches.P5571 <- silica.df[grep(paste(toMatch, collapse='|'), silica.df$Go.Cc), ]   # 72 matches
length(intersect(mem.P5571$Protein.Acs, matches.P5571$Protein.Acs)) # 41 non-PM membrane proteins

silica.df$matched_mask <- sapply( strsplit(silica.df$Go.Cc , ';', fixed = TRUE ),
                                    function ( CCs ) length( CCs ) > 0 && all( CCs %in% 'integral to membrane' ) )
silica.itm <- subset(silica.df, matched_mask)   # No matches for this subset

# Find SLCs
length(grep("SLC", silica.df$Gene))   # 1 match
slc.silica <- silica.df[grep("SLC", silica.df$Gene), ]
write.csv2(slc.amino, file="SLCs_P5571.csv")    # Save SLCs to .csv file  
  
----------------------------------------------------
######################################
# Annotate proteins from scratch via biomaRt
######################################

x <- silica.df[['Protein.Acs']]

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
                 FUN = function(X) paste(unique(X), collapse=", "))   # 270 proteins annotated
colnames(DF2) <- c("Protein.Acs", "Go.Cc")

combi <- merge(silica.df[,1:12], DF2, by="Protein.Acs")     # Combine Isobar data and annotations into one dataframe
write.csv2(combi, file="P5571.GOterms.csv")

# Read in the file created before =====================================================================================
annot <- read.csv2("~/MP_RStudio_May2014/SilicaBeads/P5571.GOterms.csv", header=TRUE, stringsAsFactors=F); dim(annot)
annot$norm.sc <- annot$Global.Spectra.Count/sum(annot$Global.Spectra.Count)   # add a column with the normalized sc

# Look for SLCs ------------------------------------------------------------
length(unique(grep("SLC", annot$Gene)))   # 1 SLCs
length(unique(grep("ABC", annot$Gene)))   # 0 ABCs
slc <- annot[grep("SLC", annot$Gene),]
abc <- annot[grep("ABC", annot$Gene),]
write.csv2(rbind(slc, abc), file="~/MP_RStudio_May2014/SilicaBeads/ABCsSLCs_P5571.csv")

# Membrane
mem <- annot[grep('membrane', annot$Go.Cc), ]; dim(mem)   # 94 proteins

# PM only
PM.only <- c('plasma membrane', 'cell surface', 'cell membrane')
matches.pm.only <- annot[grep(paste(PM.only, collapse='|'), annot$Go.Cc), ]; dim(matches.pm.only)  # 50 matches
write.csv2(matches.pm.only, file="~/MP_RStudio_May2014/SilicaBeads/PMOnly_P5571.csv")

# Search results with a vector of GOCC annotations of interest
toMatch <- c('plasma membrane', 'cell surface', 'extracellular', 'cell membrane')
matches.pm <- annot[grep(paste(toMatch, collapse='|'), annot$Go.Cc), ]   # 151 matches
length(setdiff(mem$Protein.Acs, matches.pm$Protein.Acs)) # 26 non-PM membrane proteins
non_PM_mem <- annot[grep(paste(setdiff(mem$Protein.Acs, matches.pm$Protein.Acs), collapse="|"), annot$Protein.Acs), ]
excluded <- c(matches.pm$Protein.Acs, non_PM_mem$Protein.Acs)
other <- annot[grep(paste(excluded, collapse="|"), annot$Protein.Acs, invert=T), ]
dim(matches.pm); dim(non_PM_mem); dim(other); dim(matches.pm)[1]+dim(non_PM_mem)[1]+dim(other)[1]

write.csv2(matches.pm, file="PM_P5571.csv")
write.csv2(non_PM_mem, file="Non_PM_mem_P5571.csv")
write.csv2(other, file="Other_P5571.csv")


combi$matched_mask <- sapply( strsplit(combi$Go.Cc, ';', fixed = TRUE ),
                                  function ( CCs ) length( CCs ) > 0 && all( CCs %in% 'integral to membrane' ) )
silica.itm <- subset(combi, matched_mask)   # No matches for this subset

write.csv2(matches.pm, file="PM_P5571.csv")

# Try to find out why the PM numbers differ so much for the two ways to get GO.CC annotations
check <- merge(combi[,c('Protein.Acs', "Go.Cc")], silica.df[, c("Protein.Acs", "Go.Cc")], by="Protein.Acs")
colnames(check) <- c("AC", "db2goCC", "Florian.GoCC")

length(which(check$db2goCC!=check$Florian.GoCC))    # 259 proteins have differing GoCC annotations!
diff <- check[which(check$db2goCC!=check$Florian.GoCC), ]
write.csv2(diff, file="DifferingGoCC.annotations.P5571.csv")    # The newly downloaded biomaRt annotations contain much more often the term "extracellular vesicular exosome"
