# Results for P6160. 
# Neg ctrl for P6158: Unlabelled sample, Cell Surface Protein Extraction Kit, Biotin elution, 1-D shotgun, 2 injections

setwd("~/MP_RStudio_May2014/CellSurfaceExtractionKit/P6160/"); getwd();

# Read in the file created with "P6160_results.R" containing the GOCC annotation =========================================
annot <- read.csv2("~/MP_RStudio_May2014/CellSurfaceExtractionKit/P6160//P6160.GOterms.csv", header=T, stringsAsFactors=F)
annot$norm.sc <- annot$global.spectra.count/sum(annot$global.spectra.count)   # add a column with the normalized sc

# Look for SLCs -------------------------------------
length(unique(grep("SLC", annot$Gene)))   # 1 SLCs
length(unique(grep("ABC", annot$Gene)))   # 0 ABCs

slc <- annot[grep("SLC", annot$Gene),]
abc <- annot[grep("ABC", annot$Gene),]
write.csv2(rbind(slc, abc), file="SLCsABCs_P6160.csv")  


# Divide into GOCC subgroups ----------------------------------
mem <- annot[grep('membrane', annot$Go.Cc), ]   # 9 proteins

# PM only
PM.only <- c('plasma membrane', 'cell surface', 'cell membrane')
matches.pm.only <- annot[grep(paste(PM.only, collapse='|'), annot$Go.Cc), ]; dim(matches.pm.only)  # 10 matches

# Extended PM
toMatch <- c('plasma membrane', 'cell surface', 'extracellular', 'cell membrane')
matches.pm <- annot[grep(paste(toMatch, collapse='|'), annot$Go.Cc), ]; dim(matches.pm)    # Gives 24 matches

# Non-PM membrane
length(setdiff(mem$protein.acs, matches.pm$protein.acs)) # 0 non-PM membrane proteins
non_PM_mem <- annot[grep(paste(setdiff(mem$protein.acs, matches.pm$protein.acs), collapse="|"), annot$protein.acs), ]
# Somehow, if there's not setdiff, the grep command does not work properly. Therefore, manually set non_PM_mem to NULL
non_PM_mem <- NULL

# Others
excluded <- c(matches.pm$protein.acs, non_PM_mem$protein.acs)
other <- annot[grep(paste(excluded, collapse="|"), annot$protein.acs, invert=T), ]

# Get sizes of subgroups
dim(matches.pm); dim(non_PM_mem); dim(other); dim(matches.pm)[1]+dim(non_PM_mem)[1]+dim(other)[1]

# ITMs
annot$matched_mask <- sapply( strsplit(annot$Go.Cc, ';', fixed = TRUE ),
                                   function ( CCs ) length( CCs ) > 0 && all( CCs %in% 'integral to membrane' ) )
P6160.itm <- subset(annot, matched_mask)   # No matches for this subset

write.csv2(matches.pm.only, file="PMOnly_P6160.csv")
write.csv2(matches.pm, file="PM_P6160.csv")
write.csv2(other, file="Other_P6160.csv")


# Plot sc distribution ========================================
pdf("~/MP_RStudio_May2014/CellSurfaceExtractionKit/P6160/GSC_distrib_P6160.pdf")
plot(log10(matches.pm$global.spectra.count), 
     ylab="Log10(Global Spectra Count)", 
     xlab="Proteins",
     col="red")
points(log10(non_PM_mem$global.spectra.count), col="orange")
points(log10(other$global.spectra.count))
legend("topright", legend=c("PM", "Membrane-PM", "Other"), fill=c("red", "orange", "black"))
dev.off();

# "Classic" vertical boxplot
pdf("~/MP_RStudio_May2014/CellSurfaceExtractionKit/P6160/sc_distrib_vertical_P6160.pdf")
boxplot(log10(matches.pm$global.spectra.count), #log10(non_PM_mem$global.spectra.count), 
        log10(other$global.spectra.count), 
        horizontal=F,
        ylab="Log10(Global Spectra Counts)", 
        names=c("PM", "other"),
        col=c("red", "grey"), 
        main="Global Spectra Count distribution labelled sample")
dev.off()

# Horizontal version of the boxplot
pdf("~/MP_RStudio_May2014/CellSurfaceExtractionKit/P6160/sc_distrib_horizontal_P6160.pdf")
boxplot(log10(matches.pm$global.spectra.count), #log10(non_PM_mem$global.spectra.count), 
        log10(other$global.spectra.count), 
        horizontal=T,
        ylab="Log10(Global Spectra Counts)", 
        names=c("PM", "other"),
        col=c("red", "grey"), 
        main="Global Spectra Count distribution labelled sample")
dev.off()