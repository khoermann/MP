# Making a stacked bar plot

# For subcellular localization
totals <- c(1938,1846,2036)
membrane <- c(1236,585,801)
other <- totals-membrane
vec <- rbind(membrane, other)
colnames(vec) <- c("Detergent", "Aqueous", "Control")
mat <- as.matrix(vec)

pdf("~/MP_RStudio_May2014/Triton_X-114/StackedBarplot_PhaseDistrib.pdf")
barplot(mat, col=c("Firebrick3", "darkblue"), cex.names = 1, las=1, xlim=c(0, 5.5), ylab="# of proteins")
legend("bottomleft", legend = c("Membrane", "Other"), fill = c("Firebrick3", "darkblue"))
dev.off();


# SLCs/ABCs
slc <- c(71,9,16)
abc <- c(11,2,4)
v <- rbind(slc,abc)
rownames(v) <- c("SLCs", "ABCs")
colnames(v) <- c("Detergent", "Aqueous", "Control")
m <- as.matrix(v)

# Plot with ggplot requires long data format
require(reshape2)
melt.mat <- melt(m)

require(ggplot2)
require(RColorBrewer)
lp.col <- brewer.pal(8, "Reds")
lp <- ggplot(melt.mat, aes(x=Var1, y=Var2, fill=value)) + geom_tile() + ylab("") + xlab("") 
pdf("~/MP_RStudio_May2014/Triton_X-114/Levelplot_ggplot.pdf")
lp
dev.off();