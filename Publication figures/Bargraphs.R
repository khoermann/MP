# Figure for tech paper showing total numbers of proteins
# For biotin elution, only the first two tryouts were taken as dataframes of duplicates and triplicates are not feasible

# Data assembly and reshaping =====================================================
SDS.elut <- c(881, 477)
biotin.elut <- c(440,394)
amino <- c(263,274)
silica <- c(258,129)

df <- as.matrix(cbind(SDS.elut, biotin.elut, amino, silica))
colnames(df) <- c("Sulfo-NHS-SS-Biotin SDS elution", "Sulfo-NHS-SS-Biotin Biotin elution", "Aminooxybiotin", "Silica Beads")

require(reshape2)
melt.df <- melt(df)

require(bear)
dfc <- summarySE(melt.df, measurevar="value", groupvars="X2")

# Plotting =====================================================
require(ggplot2)

# Basic color options
# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# To use for fills, add
scale_fill_manual(values=cbPalette)

# To use for line and point colors, add
scale_colour_manual(values=cbPalette)

# The plot itself
pdf("~/MP_RStudio_May2014/Publication figures/Barplot_PMTotals.pdf", paper="a4r", width=20, height=10)
ggplot(dfc, aes(x=X2, y=value))+
  theme_bw()+
  theme(axis.title.y=element_text(size=18, vjust=1.5))+
  theme(axis.text.x=element_text(size=10))+
  theme(axis.text.y=element_text(size=14))+
  ylab("PM proteins")+
  xlab("")+
  geom_bar(position=position_dodge(), stat="identity", fill="#999999", colour="black", size=.3)+
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.2, size=.5, position=position_dodge(.9))+
  ggtitle("Numbers of PM proteins obtained with different protocols")+
  theme(plot.title=element_text(size=24, vjust=1.5))+
  scale_y_continuous(breaks=1:800*100)
dev.off();


# Stacked barplot showing composition of PM fraction
totals <- c(881,477,263,274,258,129,440,394,394)
extr <- c(410,255,85,184,167,81,227,204,198)
int.PM <- c(264,113,133,23,17,14,115,103,100)
non.int.PM <- c(207,109,45,67,74,34,98,87,96)

df.stack <- rbind(extr, int.PM, non.int.PM)
colnames(df.stack) <- c(rep(c("Sulfo-NHS-SS-Biotin SDS elution", "Aminooxybiotin", "Silica Beads"), each=2), rep("Sulfo-NHS-SS-Biotin Biotin elution", 3))
rownames(df.stack) <- c("extracellular", "integral to plasma membrane", "plasma membrane")

melt.df.stack <- melt(df.stack)
bfc <- summarySE(melt.df.stack, measurevar="value", groupvars=c("X1", "X2"))

require(RColorBrewer)
col <- brewer.pal(3, "Greys")

pdf("~/MP_RStudio_May2014/Publication figures/StackedBarplot_Subcomposition_greys.pdf", paper="a4r", width=20, height=10)
ggplot(bfc, aes(x=X2, y=value, fill=X1))+
  theme_bw()+
  theme(axis.title.y=element_text(size=18, vjust=1.5))+
  theme(axis.text.x=element_text(size=10))+
  theme(axis.text.y=element_text(size=14))+
  geom_bar(stat="identity", colour="black")+
  xlab("")+
  ylab("Number of proteins")+
  ggtitle("Subsets within the purified cell surface fractions")+
  theme(plot.title=element_text(size=24, vjust=1.5))+
  scale_y_continuous(breaks=1:800*100)+
  guides(fill=guide_legend(reverse=TRUE, title=NULL))+
  scale_fill_manual(values=col, guide=guide_legend(reverse=TRUE))+
  theme(legend.position="top")+
  theme(legend.text=element_text(size=14))
dev.off();