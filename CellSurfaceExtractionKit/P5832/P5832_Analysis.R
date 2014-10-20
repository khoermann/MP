## Deeper analysis of Pierce kit data
# Modification: Biotin elution -> no FASP required

setwd('~/MP_RStudio_May2014//CellSurfaceExtractionKit/P5832/'); getwd(); 

kit.df <- read.csv2('20140325_Isobar_P5832.csv', header=TRUE, stringsAsFactors=F)
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
                 FUN = function(X) paste(unique(X), collapse=", "))   # 808 proteins annotated
colnames(DF2) <- c("protein.acs", "Go.Cc")
dim(DF2)

combi <- merge(kit.df, DF2, by="protein.acs")     # Combine Isobar data and annotations into one dataframe
write.csv2(combi, file="P5832.GOterms.csv")

# Start with file created above ============================================
annot.5832 <- read.csv2("P5832.GOterms.csv", header=T, stringsAsFactors=F)
annot.5832$norm.sc <- annot.5832$global.spectra.count/sum(annot.5832$global.spectra.count)   # add a column with the normalized sc

# Look for SLCs
length(unique(grep("SLC", annot.5832$Gene)))   # 15 SLCs
length(unique(grep("ABC", annot.5832$Gene)))   # 2 ABCs

slc.5832 <- annot.5832[grep("SLC", annot.5832$Gene),]
abc.5832 <- annot.5832[grep("ABC", annot.5832$Gene),]
write.csv2(rbind(slc.5832, abc.5832), file="SLCsABCs_P5832.csv")  

mem <- annot.5832[grep('membrane', annot.5832$Go.Cc), ]; dim(mem)  # 351 proteins

# PM only
PM.only <- c('plasma membrane', 'cell surface', 'cell membrane')
matches.pm.only <- annot.5832[grep(paste(PM.only, collapse='|'), annot.5832$Go.Cc), ]; dim(matches.pm.only)  # 213 matches
write.csv2(matches.pm.only, file="~/MP_RStudio_May2014/CellSurfaceExtractionKit//P5832/PMOnly_P5832.csv")

# Search results with a vector of GOCC annotations of interest
toMatch <- c('plasma membrane', 'cell surface', 'extracellular', 'cell membrane')
matches.pm <- annot.5832[grep(paste(toMatch, collapse='|'), annot.5832$Go.Cc), ]; dim(matches.pm)  # Gives 391 matches
length(setdiff(mem$protein.acs, matches.pm$protein.acs)) # 163 non-PM membrane proteins
non_PM_mem <- annot.5832[grep(paste(setdiff(mem$protein.acs, matches.pm$protein.acs), collapse="|"), annot.5832$protein.acs), ]
excluded <- c(matches.pm$protein.acs, non_PM_mem$protein.acs)
other <- annot.5832[grep(paste(excluded, collapse="|"), annot.5832$protein.acs, invert=T), ]
dim(matches.pm); dim(non_PM_mem); dim(other); dim(matches.pm)[1]+dim(non_PM_mem)[1]+dim(other)[1]

annot.5832$matched_mask <- sapply( strsplit(annot.5832$Go.Cc, ';', fixed = TRUE ),
                                   function ( CCs ) length( CCs ) > 0 && all( CCs %in% 'integral to membrane' ) )
P5832.itm <- subset(annot.5832, matched_mask)   # No matches for this subset

write.csv2(matches.pm, file="PM_5832.csv")
write.csv2(non_PM_mem, file="Non_PM_mem_P5832.csv")
write.csv2(other, file="Other_P5832.csv")  

=====================================================================================
  
# Plot sc distribution ========================================
pdf("~/MP_RStudio_May2014/CellSurfaceExtractionKit/P5832/GSC_distrib_P5832.pdf")
plot(log10(matches.pm$global.spectra.count), 
     ylab="Log10(Global Spectra Count)", 
     xlab="Proteins",
     col="red")
points(log10(non_PM_mem$global.spectra.count), col="orange")
points(log10(other$global.spectra.count))
legend("topright", legend=c("PM", "non-PM Mem", "Other"), fill=c("red", "orange", "black"))
dev.off();

# "Classic" vertical boxplot
pdf("~/MP_RStudio_May2014/CellSurfaceExtractionKit/P5832/sc_distrib_vertical_P5832.pdf")
boxplot(log10(matches.pm$global.spectra.count), log10(non_PM_mem$global.spectra.count), log10(other$global.spectra.count), 
        horizontal=F,
        ylab="Log10(Global Spectra Counts)", 
        names=c("PM", "Membrane-PM", "other"),
        col=c("red", "orange", "grey"), 
        main="Global Spectra Count distribution labelled sample")
dev.off()

# Horizontal version of the boxplot
pdf("~/MP_RStudio_May2014/CellSurfaceExtractionKit/P5832/sc_distrib_horizontal_P5832.pdf")
boxplot(log10(matches.pm.only$global.spectra.count), log10(non_PM_mem$global.spectra.count), log10(other$global.spectra.count), 
        horizontal=T,
        xlab="Log10(Global Spectra Counts)", 
        names=c("PM Only", "Membrane-PM", "other"), 
        col=c("red", "orange", "grey"), 
        main="Global Spectra Count distribution labelled sample")
dev.off()