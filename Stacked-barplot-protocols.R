# Making a stacked bar plot for the modified Cell Surface Protein Extraction kit tryouts

setwd('/slow/home/khoermann/MP/'); getwd();

totals <- c(2283,966,318,465,277,117,544,205)
PM <- c(968,477,256,280,151,81,258,129)
add.mem <- c(263,123,17,40,26,2,65,13)
others <- totals-PM-add.mem
vec <- rbind(totals, PM, add.mem, others)
colnames(vec) <- c('1st Biotinylation','2nd Biotinylation','1st Aminooxy', '2nd Aminooxy','1st SB (30 Mio)',
                   "2nd SB (30 Mio)", "1st SB (50 Mio)", "2nd SB (50 Mio)")
mat <- as.matrix(vec)

library(RColorBrewer)
sequential <- brewer.pal(3, "PuBu")
barplot(mat[2:4,], col=sequential, cex.names = 0.6, ylab = "# of proteins", las=1, main="Overview on results of protocol tryouts")
legend("topright", legend = c('PM', 'other membrane', 'others'), fill = sequential[1:3], title='GO annotations')

