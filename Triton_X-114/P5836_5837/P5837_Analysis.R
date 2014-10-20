# Results for P5837: 
# P5837: 30 Mio wt cells subjected to phase separation +FASP, biological replicate of P5836
# Triton X-114, aqueous phase, FASP, 1-D shotgun, 2 injections

setwd("~/MP_RStudio_May2014/Triton_X-114/P5836_5837/"); getwd();

# Read in the file created with "P5837_results.R" containing the GOCC annotation =========================================
annot <- read.csv2("~/MP_RStudio_May2014/Triton_X-114/P5836_5837/P5837.GOterms.csv", header=T, stringsAsFactors=F)
dim(annot)
annot$norm.sc <- annot$global.spectra.count/sum(annot$global.spectra.count)   # add a column with the normalized sc

# Look for SLCs -------------------------------------
length(unique(grep("SLC", annot$Gene)))   # 9 SLCs
length(unique(grep("ABC", annot$Gene)))   # 2 ABCs

slc <- annot[grep("SLC", annot$Gene),]
abc <- annot[grep("ABC", annot$Gene),]
write.csv2(rbind(slc, abc), file="SLCsABCs_P5837.csv")  


# Divide into GOCC subgroups ----------------------------------
mem <- annot[grep('membrane', annot$Go.Cc), ]   # 788 proteins
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

# Others (= non PM)
other <- annot[grep(paste(mem$protein.acs, collapse="|"), annot$protein.acs, invert=T), ]
dim(mem); dim(matches.pm); dim(non_PM_mem); dim(other); 

# ITMs
annot$matched_mask <- sapply( strsplit(annot$Go.Cc, ';', fixed = TRUE ),
                                   function ( CCs ) length( CCs ) > 0 && all( CCs %in% 'integral to membrane' ) )
P5837.itm <- subset(annot, matched_mask)   # No matches for this subset
dim(P5837.itm)

write.csv2(matches.pm.only, file="PMOnly_P5837.csv")
write.csv2(matches.pm, file="PM_P5837.csv")
write.csv2(non_PM_mem, file="Non_PM_mem_P5837.csv")
write.csv2(other, file="Other_P5837.csv")
write.csv2(mem, file="Membrane_P5837.csv")


# Plot sc distribution ========================================
pdf("~/MP_RStudio_May2014/Triton_X-114/P5834_5835/GSC_distrib_P5837.pdf")
plot(log10(mem$global.spectra.count), 
     ylab="Log10(Global Spectra Count)", 
     xlab="Proteins",
     col="red")
points(log10(other$global.spectra.count))
legend("topright", legend=c("Membrane", "Non-membrane"), fill=c("red", "black"))
dev.off();

# "Classic" vertical boxplot
pdf("~/MP_RStudio_May2014/Triton_X-114/P5834_5835//sc_distrib_vertical_P5837.pdf")
boxplot(log10(mem$global.spectra.count), log10(other$global.spectra.count), 
        horizontal=F,
        ylab="Log10(Global Spectra Counts)", 
        names=c("Membrane", "Non-membrane"),
        col=c("red", "grey"), 
        main="Global Spectra Count distribution P5837")
dev.off()

# Horizontal version of the boxplot
pdf("~/MP_RStudio_May2014/Triton_X-114/P5834_5835//sc_distrib_horizontal_P5837.pdf")
boxplot(log10(mem$global.spectra.count), log10(other$global.spectra.count), 
        horizontal=T,
        ylab="Log10(Global Spectra Counts)", 
        names=c("Membrane", "Non-membrane"),
        col=c("red", "grey"), 
        main="Global Spectra Count distribution P5837")
dev.off()