## 2nd tryout Pierce kit

setwd('~/MP_RStudio_May2014/CellSurfaceExtractionKit/P5826/'); getwd(); 

kit.df <- read.csv2('20140325_Isobar_P5826.csv', header=TRUE, stringsAsFactors=F)   # 966 proteins
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
                 FUN = function(X) paste(unique(X), collapse=", "))   # 949 proteins annotated
colnames(DF2) <- c("protein.acs", "Go.Cc")

combi <- merge(kit.df, DF2, by="protein.acs")     # Combine Isobar data and annotations into one dataframe
write.csv2(combi, file="P5826.GOterms.csv")
  
annot.5826 <- read.csv2("~/MP_RStudio_May2014/CellSurfaceExtractionKit/P5826/P5826.GOterms.csv", header=T, stringsAsFactors=F)
annot.5826$norm.sc <- annot.5826$global.spectra.count/sum(annot.5826$global.spectra.count)   # add a column with the normalized sc

# Look for SLCs
length(unique(grep("SLC", annot.5826$Gene)))   # 17 SLCs
length(unique(grep("ABC", annot.5826$Gene)))   # 1 ABCs

slc.5826 <- annot.5826[grep("SLC", annot.5826$Gene),]
abc.5826 <- annot.5826[grep("ABC", annot.5826$Gene),]
write.csv2(rbind(slc.5826, abc.5826), file="SLCsABCs_P5826.csv")  

mem <- annot.5826[grep('membrane', annot.5826$Go.Cc), ]   # 393 proteins

# PM only
PM.only <- c('plasma membrane', 'cell surface', 'cell membrane')
matches.pm.only <- annot.5826[grep(paste(PM.only, collapse='|'), annot.5826$Go.Cc), ]; dim(matches.pm.only)  # 222 matches
write.csv2(matches.pm.only, file="~/MP_RStudio_May2014/CellSurfaceExtractionKit//P5826/PMOnly_P5826.csv")

# Search results with a vector of GOCC annotations of interest
toMatch <- c('plasma membrane', 'cell surface', 'extracellular', 'cell membrane')
matches.pm <- annot.5826[grep(paste(toMatch, collapse='|'), annot.5826$Go.Cc), ]    # Gives 477 matches
length(setdiff(mem$protein.acs, matches.pm$protein.acs)) # 123 non-PM membrane proteins
non_PM_mem <- annot.5826[grep(paste(setdiff(mem$protein.acs, matches.pm$protein.acs), collapse="|"), annot.5826$protein.acs), ]
excluded <- c(matches.pm$protein.acs, non_PM_mem$protein.acs)
other <- annot.5826[grep(paste(excluded, collapse="|"), annot.5826$protein.acs, invert=T), ]
dim(matches.pm); dim(non_PM_mem); dim(other); dim(matches.pm)[1]+dim(non_PM_mem)[1]+dim(other)[1]

annot.5826$matched_mask <- sapply( strsplit(annot.5826$Go.Cc, ';', fixed = TRUE ),
                              function ( CCs ) length( CCs ) > 0 && all( CCs %in% 'integral to membrane' ) )
P5826.itm <- subset(annot.5826, matched_mask)   # No matches for this subset

write.csv2(matches.pm, file="PM_5826.csv")
write.csv2(non_PM_mem, file="Non_PM_mem_P5826.csv")
write.csv2(other, file="Other_P5826.csv")

# Plot sc distribution ========================================
pdf("~/MP_RStudio_May2014/CellSurfaceExtractionKit/P5826/GSC_distrib.pdf")
plot(log10(matches.pm$global.spectra.count), 
     ylab="Log10(Global Spectra Count)", 
     xlab="Proteins",
     col="red")
points(log10(non_PM_mem$global.spectra.count), col="orange")
points(log10(other$global.spectra.count))
legend("topright", legend=c("PM", "non-PM Mem", "Other"), fill=c("red", "orange", "black"))
dev.off();

# "Classic" vertical boxplot
pdf("~/MP_RStudio_May2014/CellSurfaceExtractionKit/P5826/sc_distrib_vertical.pdf")
boxplot(log10(matches.pm$global.spectra.count), log10(non_PM_mem$global.spectra.count), log10(other$global.spectra.count), 
        horizontal=F,
        ylab="Log10(Global Spectra Counts)", 
        names=c("PM", "Mem non-PM", "other"), 
        col=c("red", "orange", "grey"), 
        main="Global Spectra Count distribution labelled sample")
dev.off()

# Horizontal version of the boxplot
pdf("~/MP_RStudio_May2014/CellSurfaceExtractionKit/P5826/sc_distrib_horizontal.pdf")
boxplot(log10(matches.pm.only$global.spectra.count), log10(non_PM_mem$global.spectra.count), log10(other$global.spectra.count), 
        horizontal=T,
        xlab="Log10(Global Spectra Counts)", 
        names=c("PM Only", "Mem non-PM", "other"), 
        col=c("red", "orange", "grey"), 
        main="Global Spectra Count distribution labelled sample")
dev.off()

###############################################################
# Alternative: Search with the Accession codes from Swissprot
sp.pm.acc <- sp.pm$V1
overlap <- intersect(sp.pm.acc, kit.df$accession)   # ONLY 103
kit.sp <- kit.df[grep(paste(overlap, collapse='|'), kit.df$accession), ]

sp.pmcsext <- read.table('Swissprot_PMCSEXT_human_ACC_20140328.txt', stringsAsFactors=F)    # Contains 6418 Protein ACs (28.03.2014)
sp.pmcsext.acc <- sp.pmcsext$V1     # ONLY 108
overlap.pmcsext <- intersect(sp.pmcsext.acc, kit.df$accession)
# Identify 'contaminations' by other membrane-bound organelles
length(grep(paste(c('nuclear', 'nucleus'), collapse='|'), matches.not.pm$go_name))   # absolute
length(grep(paste(c('nuclear', 'nucleus'), collapse='|'), matches.not.pm$go_name))/dim(kit.df)[1]    # relative
post.nuclear <- matches.not.pm[grep(paste(c('nuclear', 'nucleus'), collapse='|'), matches.not.pm$go_name, invert=TRUE), ]

length(grep(paste(c('cytosol', 'cytoplasm'), collapse='|'), post.nuclear$go_name))   # absolute
length(grep(paste(c('cytosol', 'cytoplasm'), collapse='|'), post.nuclear$go_name))/dim(kit.df)[1]    # relative
post.nuclear.cyto <- post.nuclear[grep(paste(c('cytosol', 'cytoplasm'), collapse='|'), post.nuclear$go_name, invert=TRUE), ]

length(grep("endoplasmic reticulum", post.nuclear.cyto$go_name))   # absolute
length(grep("endoplasmic reticulum", post.nuclear.cyto$go_name))/dim(kit.df)[1]    # relative

length(grep("Golgi", post.nuclear.cyto$go_name))   # absolute
length(grep("Golgi", post.nuclear.cyto$go_name))/dim(kit.df)[1]    # relative

length(grep(paste(c("mitochondria","mitochondrial"), collapse="|"), post.nuclear.cyto$go_name))   # absolute
length(grep(paste(c("mitochondria","mitochondrial"), collapse="|"), post.nuclear.cyto$go_name))/dim(kit.df)[1]   # absolute
