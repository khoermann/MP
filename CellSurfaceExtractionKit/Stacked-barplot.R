# Making a stacked bar plot for the modified Cell Surface Protein Extraction kit tryouts

setwd('~/MP_RStudio_May2014//CellSurfaceExtractionKit/'); getwd();

totals <- c(795,966,821,574,648)
PM.ext <- c(391,477,440,281,347)
PM.only <- c(165,222,213,136,158)
ext <- PM.ext-PM.only
add.mem <- c(163,123,107,72,77)
others <- totals-PM.ext-add.mem
vec <- rbind(PM.ext, add.mem, others)
colnames(vec) <- c('deglycosylated','classic','Biotin elution','+TX-114 detergent', '+TX-114 acqueous')
mat <- as.matrix(vec)

library(RColorBrewer)
sequential <- brewer.pal(4, "Blues")
pdf("~/MP_RStudio_May2014/CellSurfaceExtractionKit//pdf/StackedBarplot.pdf")
barplot(mat, col=sequential, cex.names = 0.6, ylab = "# of proteins", ylim=c(0,1000), xlim=c(0,8), las=1)
legend("bottomright", legend = c('PM', 'PM extended', 'other membrane', 'others'), fill = sequential[1:4], title='GO annotations')
dev.off();

# "Small version" just containing the first three columns =========================================
barplot(mat[,1:3], col=sequential, cex.names = 0.8, ylab = "# of proteins", ylim=c(0,1000), xlim=c(0,10), las=1)
legend("bottom", legend = c('PM', 'Membrane-PM', 'others'), fill = sequential[1:4], title='GO annotations')



# Include 'negative controls'
totals.neg <- c(795,546,966,606,574,544,648,414,821,341)
PM.neg <- c(391,266,482,287,281,349,347,205,441,169)
add.mem.neg <- c(167,94,164,85,99,65,90,46,129,45)
others.neg <- totals.neg-PM.neg-add.mem.neg
vec.neg <- rbind(totals.neg, PM.neg, add.mem.neg, others.neg)
colnames(vec.neg) <- c('deglycosylated', 'neg ctrl','classic','neg ctrl','+TX-114 detergent','neg ctrl', '+TX-114 acqueous','neg ctrl','Biotin elution','neg ctrl')
mat.neg <- as.matrix(vec.neg)

library(RColorBrewer)
# This part needs editing
barplot(mat.neg[2:4,], col=sequential, cex.names = 0.6, ylab = "# of proteins", ylim=c(0,1000), xlim=c(0,16), las=1)
legend("bottomright", legend = c('PM', 'other membrane', 'others'), fill = sequential[1:3], title='GO annotations')
