# Making a stacked bar plot

# For the Cell Surface Extraction kit
PM <- c(31, 50)
Nucleus <- c(22,22)
Cytoplasma <- c(20,10)
Mitochondria <- c(3,6)
ER <- c(3,3)
Golgi <- c(2,1)
others <- 100-PM-ER-Golgi-Nucleus-Mitochondria-Cytoplasma
vec <- rbind(PM, Nucleus, Cytoplasma, Mitochondria, ER, Golgi, others)
colnames(vec) <- c("1st tryout", "2nd tryout")
mat <- as.matrix(vec)

library(RColorBrewer)
sequential <- brewer.pal(7, "PuBu")
barplot(mat, col=sequential, cex.names = 1, ylab = "%", xlim=c(0,3.5), las=1)
legend("bottomright", legend = c('others', 'Golgi', 'ER', 'Mitochondria', 'Cytoplasma', 'Nucleus', 'PM'), fill = sequential[7:1], title='Subcellular localization')


# Aminooxy approach
PM <- c(63,61)
Nucleus <- c(8,20)
Cytoplasma <- c(16,9)
Mitochondria <- c(3,3)
ER <- c(0.9,2)
Golgi <- c(0.6,0.2)
others <- 100-PM-ER-Golgi-Nucleus-Mitochondria-Cytoplasma
vec <- rbind(PM, Nucleus, Cytoplasma, Mitochondria, ER, Golgi, others)
colnames(vec) <- c("1st tryout", "2nd tryout")
mat <- as.matrix(vec)

library(RColorBrewer)
sequential <- brewer.pal(7, "PuRd")
barplot(mat, col=sequential, cex.names = 1, ylab = "%", xlim=c(0,3.5), las=1)
legend("bottomright", legend = c('others', 'Golgi', 'ER', 'Mitochondria', 'Cytoplasma', 'Nucleus', 'PM'), fill = sequential[7:1], title='Subcellular localization')

