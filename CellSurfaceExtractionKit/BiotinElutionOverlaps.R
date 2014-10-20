# Comparing the samples obtained with Biotin elution
# P5832 + P6157

# Totals
total.P5832 <- read.csv2("~/MP_RStudio_May2014/CellSurfaceExtractionKit/P5832/20140325_Isobar_P5832.csv", header=T, stringsAsFactors=F)
total.P6157 <- read.csv2("~/MP_RStudio_May2014/CellSurfaceExtractionKit/P6157/Isobar_P6157.csv", header=T, stringsAsFactors=F)

P5832.ac <- gsub("[,-].*",'',total.P5832$protein.acs)
P6157.ac <- gsub("[,-].*",'',total.P6157$protein.acs)

length(P5832.ac); length(P6157.ac)
length(intersect(P5832.ac, P6157.ac))

# Graphical output
library(VennDiagram);
library(extrafont); 
grid.newpage();
venn.plot <- draw.pairwise.venn(821,729,589, c('First','Second'), cex=rep(2.5,3),
                                col=c('plum1', 'plum3'), fill=c('plum1', 'plum3'), 
                                cat.pos=c(315,135), cat.dist=c(0.05,0.08), 
                                cat.fontfamily=rep('Franklin Gothic Book',2), cat.cex=rep(2.5,2), 
                                cat.col=c('plum1', 'plum3'),
                                scaled=TRUE, margin=c(0.08, 0.08), fontfamily=rep('Franklin Gothic Book',3));



# PM ext ================================
pm.P5832 <- read.csv2("~/MP_RStudio_May2014/CellSurfaceExtractionKit/P5832/PM_5832.csv", header=T, stringsAsFactors=F)
pm.P6157 <- read.csv2("~/MP_RStudio_May2014/CellSurfaceExtractionKit/P6157/PM_6157.csv", header=T, stringsAsFactors=F)

length(pm.P5832$protein.acs); length(pm.P6157$protein.acs)
length(intersect(pm.P5832$protein.acs, pm.P6157$protein.acs))

grid.newpage();
venn.plot <- draw.pairwise.venn(440,394,340, c('First','Second'), cex=rep(2.5,3),
                                col=c('plum1', 'plum3'), fill=c('plum1', 'plum3'), 
                                cat.pos=c(315,135), cat.dist=c(0.05,0.08), 
                                cat.fontfamily=rep('Franklin Gothic Book',2), cat.cex=rep(2.5,2), 
                                cat.col=c('plum1', 'plum3'),
                                scaled=TRUE, margin=c(0.08, 0.08), fontfamily=rep('Franklin Gothic Book',3));
