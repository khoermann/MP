# Script to intersect the proteins detected with different phases of Triton x-114 phase separation protocol
setwd("~/MP_RStudio_May2014/Triton_X-114/"); getwd();

# Read in the main data ==========================================================
det <- read.csv2("~/MP_RStudio_May2014/Triton_X-114//P5834_5835/20140325_Isobar_P5834.P5835.csv", header=T, stringsAsFactors=F)
aqu <- read.csv2("~/MP_RStudio_May2014/Triton_X-114//P5836_5837/20140325_Isobar_P5836.P5837.csv", header=T, stringsAsFactors=F)
ctrl <- read.csv2("~/MP_RStudio_May2014/Triton_X-114//P5838_5839/20140325_Isobar_P5838.P5839.csv", header=T, stringsAsFactors=F)
dim(det); dim(aqu); dim(ctrl);

require(xlsx)
det.first <- read.xlsx2("~/MP_RStudio_May2014/Triton_X-114//P5574/20131202.TritonX114.P5574.results.xlsx", sheetIndex=1, stringsAsFactors=F)

det.ac <- gsub("[,-].*",'',det$protein.acs)
aqu.ac <- gsub("[,-].*",'',aqu$protein.acs)
ctrl.ac <- gsub("[,-].*",'',ctrl$protein.acs)

det.first.ac <- gsub("[,-].*",'',det.first$Protein.Acs)

# Determine the overlap
Reduce(intersect, list(det.ac, aqu.ac, ctrl.ac))
length(intersect(det.ac, det.first.ac))
length(intersect(det.ac, aqu.ac))

# Display it graphically
library(gplots)
pdf("~/MP_RStudio_May2014/Triton_X-114/pdf/Overlaps_allProteins_det.aqu.ctrl.pdf")
venn(list(det=det.ac,
          aqu=aqu.ac,
          ctrl=ctrl.ac))
dev.off();

# Control this result again by "summing" det and aqu phase and compare it to ctrl 
extr <- unique(c(det.ac, aqu.ac))
length(intersect(extr, unique(ctrl.ac)))
venn(list(extr=extr,
          ctrl=unique(ctrl.ac)))


# Read in the membrane subsets
Mem.det <- read.csv2("~/MP_RStudio_May2014/Triton_X-114//P5834_5835/Mem_P5834.P5835.csv", header=T, stringsAsFactors=F)
Mem.aqu <- read.csv2("~/MP_RStudio_May2014/Triton_X-114//P5836_5837/Mem_P5836.P5837.csv", header=T, stringsAsFactors=F)
Mem.ctrl <- read.csv2("Mem_P5838.P5839.csv", header=T, stringsAsFactors=F)

library(gplots)
pdf("Overlaps.pdf")
venn(list(detergent=Mem.det$protein.acs,
          aqueous=Mem.aqu$protein.acs,
          control=Mem.ctrl$protein.acs))
dev.off();

# Max number of membrane proteins obtained from two phases
length(intersect(Mem.det$protein.acs, Mem.aqu$protein.acs))

# Show reproducibility: Compare to det phase first tryout P5574
Mem.det.first <- read.csv2("~/MP_RStudio_May2014/Triton_X-114/P5574//Mem_P5574.csv", header=T, stringsAsFactors=F)
length(intersect(Mem.det.first$Protein.Acs, Mem.det$protein.acs))

# If the overlap between two sets is shown, this works better --------------------------------
# Graphical output
library(VennDiagram);
library(extrafont); 
grid.newpage();
venn.plot <- draw.pairwise.venn(1938,1846,721, c('Detergent','Aqueous'), cex=rep(3.5,3),
                                col=c('red3', 'royalblue4'), fill=c('red3', 'royalblue4'), 
                                cat.pos=c(315,135), cat.dist=c(0.08,0.07), 
                                cat.fontfamily=rep('Franklin Gothic Book',2), cat.cex=rep(2.5,2), 
                                cat.col=c('red3', 'royalblue4'),
                                scaled=TRUE, margin=c(0.08, 0.08), fontfamily=rep('Franklin Gothic Book',3));


