## Checking the overlaps for the protocols

# Read in the PM proteins identified with Pierce Kit and modifications of it
PM.P5826 <- read.csv2('PM_P5826.csv', header=TRUE, stringsAsFactors=F)
head(PM.P5826)
str(PM.P5826)

PM.5824 <- read.csv2('PM_P5824.csv', header=TRUE, stringsAsFactors=F)
PM.5828 <- read.csv2('PM_P5828.csv', header=TRUE, stringsAsFactors=F)
PM.5830 <- read.csv2('PM_P5830.csv', header=TRUE, stringsAsFactors=F)
PM.5832 <- read.csv2('PM_P5832.csv', header=TRUE, stringsAsFactors=F)

library(gplots)
test1 <- PM.5824$accession
test2 <- PM.P5826$accession
test3 <- PM.5828$accession
test4 <- PM.5830$accession
test5 <- PM.5832$accession

venn(list(deglyco=test1,
          classic=test2,
          detergent=test3,
          acqueous=test4, 
          Biotin=test5))


# Comparisons within the other, 'established' protocols
# Cell Surface Protein Extraction Kit
P5134 <- read.csv2("~/MP_RStudio_May2014/CellSurfaceExtractionKit/P5134-5136/PM_P5134.csv", 
                        header=T, stringsAsFactors=F); dim(P5134);   # 881 PM proteins
P5826 <- read.csv2("~/MP_RStudio_May2014/CellSurfaceExtractionKit/P5826/PM_5826.csv", 
                   header=T, stringsAsFactors=F); dim(P5826);   # 477 PM proteins
length(intersect(P5134$protein.acs, P5826$protein.acs))    # Overlap: 398 proteins

      
# Aminooxybiotin
P5455 <- read.csv2("~/MP_RStudio_May2014/Aminooxybiotin/P5396_P5455/P5455/PM_P5455.csv", 
                        header=T, stringsAsFactors=F); dim(P5455);   # 263 PM proteins
P5767 <- read.csv2("~/MP_RStudio_May2014/Aminooxybiotin//P5767_P5768/P5767/PM_P5767.csv", 
                        header=T, stringsAsFactors=F) # 274 PM proteins
length(intersect(P5455$protein.acs, P5767$protein.acs))   # Overlap: 107
P5396 <- read.csv2("~/MP_RStudio_May2014/Aminooxybiotin/P5396_P5455/P5396/PM_P5396.csv", 
                   header=T, stringsAsFactors=F); dim(P5396);   # 37 PM proteins
P5768 <- read.csv2("~/MP_RStudio_May2014/Aminooxybiotin//P5767_P5768/P5768/PM_P5768.csv", 
                   header=T, stringsAsFactors=F) # 255 PM proteins
length(intersect(P5455$protein.acs, P5396$protein.acs))   # Overlap: 34
length(intersect(P5767$protein.acs, P5768$protein.acs))   # Overlap: 210
length(intersect(P5455$protein.acs, P5768$protein.acs))   # Overlap: 210

# Silica Beads
# 30 Mio cells
P5571 <- read.csv2("~/MP_RStudio_May2014/SilicaBeads/PM_P5571.csv", header=T, stringsAsFactors=F) # 151 PM proteins
P5744 <- read.csv2("~/MP_RStudio_May2014/SilicaBeads/PM_P5744.csv", header=T, stringsAsFactors=F) # 81 PM proteins
length(intersect(P5571$Protein.Acs, P5744$protein.acs))   # Overlap: 39 

# 50 Mio cells
P5572 <- read.csv2("~/MP_RStudio_May2014/SilicaBeads/PM_P5572.csv", header=T, stringsAsFactors=F) # 258 proteins
P5745 <- read.csv2("~/MP_RStudio_May2014/SilicaBeads/PM_P5745.csv", header=T, stringsAsFactors=F) # 129 proteins
length(intersect(P5572$Protein.Acs, P5745$protein.acs)) # Overlap: 87

# Triton Phase separation
P5834_5835 <- read.csv2("Mem_P5834.P5835.csv", header=T, stringsAsFactors=F); dim(P5834_5835);  # 1180 proteins
P5574 <- read.csv2("Mem_P5574.csv", header=T, stringsAsFactors=F); dim(P5574);  # 733 proteins
length(intersect(P5834_5835$protein.acs, P5574$Protein.Acs))  # Overlap: 697

# Modified Cell Surface Protein Extraction Kit +TX-114 phase separation
PM_5828 <- read.csv2("/slow/home/khoermann/MP/CellSurfaceExtractionKit/P5828/PM_5828.csv") # detergent
PM_5830 <- read.csv2("/slow/home/khoermann/MP/CellSurfaceExtractionKit/P5830/PM_5830.csv") # aqueous
length(intersect(PM_5828$protein.acs, PM_5830$protein.acs))   # Overlap between aqueous and detergent is 182

slc_5828 <- read.csv2("/slow/home/khoermann/MP/CellSurfaceExtractionKit/P5828/SLCsABCs_P5828.csv") # detergent
slc_5830 <- read.csv2("/slow/home/khoermann/MP/CellSurfaceExtractionKit/P5830/SLCsABCs_P5830.csv") # aqueous 
length(intersect(slc_5828$protein.acs, slc_5830$protein.acs))

library(gplots)
venn(list(slc_detergent=slc_5828$protein.acs,
          slc_aqueous=slc_5830$protein.acs))
----------------------------------------------

# Graphical output
library(VennDiagram);
library(extrafont); 
grid.newpage();
venn.plot <- draw.pairwise.venn(730,1236,690, c('First','Second'), cex=rep(2.5,3),
                                col=c('yellowgreen', 'olivedrab3'), fill=c('yellowgreen', 'olivedrab3'), 
                                cat.pos=c(315,135), cat.dist=c(0.05,0.08), 
                                cat.fontfamily=rep('Franklin Gothic Book',2), cat.cex=rep(2.5,2), 
                                cat.col=c('yellowgreen', 'olivedrab3'),
                                scaled=T, margin=c(0.08, 0.08), fontfamily=rep('Franklin Gothic Book',3));
