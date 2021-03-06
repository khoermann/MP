# Results for P6287: 
# P6287: 40 Mio wt cells labeled, Cell Surface Extraction Kit, washing with regular kit washing buffer, Biotin elution
# Cell Surface Protein Extraction Kit, Biotin elution, 1-D shotgun, 2 injections

setwd("~/MP_RStudio_May2014/CellSurfaceExtractionKit/P6287/"); getwd();

# Read in the file created with "P6287_results.R" containing the GOCC annotation =========================================
annot <- read.csv2("~/MP_RStudio_May2014/CellSurfaceExtractionKit/P6287//P6287.GOterms.csv", header=T, stringsAsFactors=F)
dim(annot)
annot$norm.sc <- annot$global.spectra.count/sum(annot$global.spectra.count)   # add a column with the normalized sc

# Look for SLCs -------------------------------------
length(unique(grep("SLC", annot$Gene)))   # 11 SLCs
length(unique(grep("ABC", annot$Gene)))   # 0 ABCs

slc <- annot[grep("SLC", annot$Gene),]
abc <- annot[grep("ABC", annot$Gene),]
write.csv2(rbind(slc, abc), file="SLCsABCs_P6287.csv")  


# Divide into GOCC subgroups ----------------------------------
mem <- annot[grep('membrane', annot$Go.Cc), ]   # 299 proteins
dim(mem)

# PM only
PM.only <- c('plasma membrane', 'cell surface', 'cell membrane')
matches.pm.only <- annot[grep(paste(PM.only, collapse='|'), annot$Go.Cc), ]; dim(matches.pm.only)  # 173 matches

# Extended PM
toMatch <- c('plasma membrane', 'cell surface', 'extracellular', 'cell membrane')
matches.pm <- annot[grep(paste(toMatch, collapse='|'), annot$Go.Cc), ]; dim(matches.pm)   # Gives 380 matches

# Non-PM membrane
length(setdiff(mem$protein.acs, matches.pm$protein.acs)) # 85 non-PM membrane proteins
non_PM_mem <- annot[grep(paste(setdiff(mem$protein.acs, matches.pm$protein.acs), collapse="|"), annot$protein.acs), ]

# Others
excluded <- c(matches.pm$protein.acs, non_PM_mem$protein.acs)
other <- annot[grep(paste(excluded, collapse="|"), annot$protein.acs, invert=T), ]
dim(matches.pm); dim(non_PM_mem); dim(other); dim(matches.pm)[1]+dim(non_PM_mem)[1]+dim(other)[1]

# ITMs
annot$matched_mask <- sapply( strsplit(annot$Go.Cc, ';', fixed = TRUE ),
                                   function ( CCs ) length( CCs ) > 0 && all( CCs %in% 'integral to membrane' ) )
P6287.itm <- subset(annot, matched_mask)   # No matches for this subset
dim(P6287.itm)

write.csv2(matches.pm.only, file="PMOnly_P6287.csv")
write.csv2(matches.pm, file="PM_P6287.csv")
write.csv2(non_PM_mem, file="Non_PM_mem_P6287.csv")
write.csv2(other, file="Other_P6287.csv")


# Plot sc distribution ========================================
pdf("~/MP_RStudio_May2014/CellSurfaceExtractionKit/P6287/GSC_distrib_P6287.pdf")
plot(log10(matches.pm$global.spectra.count), 
     ylab="Log10(Global Spectra Count)", 
     xlab="Proteins",
     col="red")
points(log10(non_PM_mem$global.spectra.count), col="orange")
points(log10(other$global.spectra.count))
legend("topright", legend=c("PM", "non-PM Mem", "Other"), fill=c("red", "orange", "black"))
dev.off();

# "Classic" vertical boxplot
pdf("~/MP_RStudio_May2014/CellSurfaceExtractionKit/P6287/sc_distrib_vertical_P6287.pdf")
boxplot(log10(matches.pm$global.spectra.count), log10(non_PM_mem$global.spectra.count), log10(other$global.spectra.count), 
        horizontal=F,
        ylab="Log10(Global Spectra Counts)", 
        names=c("PM", "Membrane-PM", "other"),
        col=c("red", "orange", "grey"), 
        main="Global Spectra Count distribution labelled sample")
dev.off()

# Horizontal version of the boxplot with PM only
pdf("~/MP_RStudio_May2014/CellSurfaceExtractionKit/P6287/sc_distrib_horizontal_P6287.pdf")
boxplot(log10(matches.pm.only$global.spectra.count), log10(non_PM_mem$global.spectra.count), log10(other$global.spectra.count), 
        horizontal=T,
        xlab="Log10(Global Spectra Counts)", 
        names=c("PM only", "Membrane-PM", "other"), 
        col=c("red", "orange", "grey"), 
        main="Global Spectra Count distribution labelled sample")
dev.off()