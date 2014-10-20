## Checking the overlaps for the three tryouts of Biotin elution

# Read in the PM proteins identified with Pierce Kit and modifications of it
PM.5832 <- read.csv2('~/MP_RStudio_May2014/CellSurfaceExtractionKit/P5832/PM_5832.csv', header=TRUE, stringsAsFactors=F)
head(PM.5832)
str(PM.5832)

PM.6157 <- read.csv2('~/MP_RStudio_May2014/CellSurfaceExtractionKit/P6157/PM_6157.csv', header=TRUE, stringsAsFactors=F)
PM.6287 <- read.csv2('~/MP_RStudio_May2014/CellSurfaceExtractionKit/P6287//PM_P6287.csv', header=TRUE, stringsAsFactors=F)
PM.6284 <- read.csv2('~/MP_RStudio_May2014/CellSurfaceExtractionKit/P6284//PM_P6284.csv', header=TRUE, stringsAsFactors=F)

# Preliminary Venn Diagrams
library(gplots)
test1 <- PM.5832$protein.acs
test2 <- PM.6157$protein.acs
test3 <- PM.6287$protein.acs
test4 <- PM.6284$protein.acs

pdf("~/MP_RStudio_May2014/CellSurfaceExtractionKit/Intersect_BiotinTryouts.pdf")
venn(list(first=test1,
          second=test2,
          third=test3,
          Urea=test4))
dev.off();

pdf("~/MP_RStudio_May2014/CellSurfaceExtractionKit/Intersect_BiotinTryouts_withoutUrea.pdf")
venn(list(first=test1,
          second=test2,
          third=test3))
dev.off();

# Get a list of the shared proteins
shared <- Reduce(intersect, list(test1, test2, test3))

# Extract spectral count info for the shared ones
shared.5832 <- PM.5832[grep(paste(shared, collapse="|"), PM.5832$protein.acs), c("protein.acs", "ID", "Description", "Gene", "global.spectra.count", "norm.sc")]
shared.6157 <- PM.6157[grep(paste(shared, collapse="|"), PM.6157$protein.acs), c("protein.acs", "ID", "Description", "Gene", "global.spectra.count", "norm.sc")]
shared.6287 <- PM.6287[grep(paste(shared, collapse="|"), PM.6287$protein.acs), c("protein.acs", "ID", "Description", "Gene", "global.spectra.count", "norm.sc")]

par(mfrow=c(1,3))
plot(shared.5832$global.spectra.count, ylim=c(0,120))
plot(shared.6157$global.spectra.count, ylim=c(0,120))
plot(shared.6287$global.spectra.count, ylim=c(0,120))
     
# Create one dataframe holding all 277 shared proteins and their gsc data
m1 <- merge(shared.5832, shared.6157[, c("protein.acs", "global.spectra.count", "norm.sc")], by="protein.acs")
m2 <- merge(m1, shared.6287[, c("protein.acs", "global.spectra.count", "norm.sc")], by="protein.acs")
colnames(m2) <- c("AC", "ID", "Description", "Gene", "gsc.5832", "norm.sc.5832", "gsc.6157","norm.sc.6157", "gsc.6287", "norm.sc.6287")

# Make one df with gsc and one with norm.sc to make statistics easier
m3 <- m2[, c(1:5, 7,9)]
m4 <- m2[, c(1:4,6,8,10)]

# Remove those where scs are 0
rm.m3 <- subset(m3, subset=(gsc.5832!=0 & gsc.6157!=0 & gsc.6287!=0))
rm.m4 <- subset(m4, subset=(norm.sc.5832!=0 & norm.sc.6157!=0 & norm.sc.6287!=0))

# Reshape to be able to apply summarySE to get statistics over all three biological replicates =======================
require(reshape2)
mm3 <- melt.data.frame(rm.m3)
mm3$ID <- as.factor(mm3$ID)

mm4 <- melt.data.frame(rm.m4)
mm4$ID <- as.factor(mm4$ID)

require(bear)
summ3 <- summarySE(mm3, measurevar="value", groupvars="ID")
summ4 <- summarySE(mm4, measurevar="value", groupvars="ID")

# Plotting part =============================================
require(ggplot2)

# Plotting sds doesn't really give nice, easily interpretable graphs
ggplot(summ3, aes(x=ID, y=sd))+
  geom_point()

ggplot(summ4, aes(x=ID, y=sd))+
  geom_point()

# Rather stay with simple correlation graphs +++++++++++++++++++++++++++++++++++++++++++++
# Use log10 transform on both axes to make readibility easier
plot1 <- ggplot(rm.m3, aes(x=gsc.5832, y=gsc.6157))+
  theme_bw()+
  scale_y_log10()+
  scale_x_log10()+
  ylab("global spectra count Rep. 2")+
  xlab("global spectra count Rep. 1")+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  theme(axis.text.x=element_text(size=14))+
  theme(axis.text.y=element_text(size=14))+
  geom_point(size=3)+
  geom_abline(intercept=0, slope=1, colour="red", size=1)+
  geom_text(aes(x=15, y=100, label=paste("r^2 ==", cor(rm.m3$gsc.5832, rm.m3$gsc.6157), sep="")), size=5, parse=TRUE)+
  theme(legend.position="none")

plot2 <- ggplot(rm.m3, aes(x=gsc.5832, y=gsc.6287))+
  theme_bw()+
  scale_y_log10()+
  scale_x_log10()+
  ylab("global spectra count Rep. 3")+
  xlab("global spectra count Rep. 1")+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  theme(axis.text.x=element_text(size=14))+
  theme(axis.text.y=element_text(size=14))+
  geom_point(size=3)+
  geom_abline(intercept=0, slope=1, colour="red", size=1)+
  geom_text(aes(x=15, y=100, label=paste("r^2 ==", cor(rm.m3$gsc.5832, rm.m3$gsc.6287), sep="")), size=5, parse=TRUE)+
  theme(legend.position="none")

plot3 <- ggplot(rm.m3, aes(x=gsc.6157, y=gsc.6287))+
  theme_bw()+
  scale_y_log10()+
  scale_x_log10()+
  ylab("global spectra count Rep. 3")+
  xlab("global spectra count Rep. 2")+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  theme(axis.text.x=element_text(size=14))+
  theme(axis.text.y=element_text(size=14))+
  geom_point(size=3)+
  geom_abline(intercept=0, slope=1, colour="red", size=1)+
  geom_text(aes(x=15, y=100, label=paste("r^2 ==", cor(rm.m3$gsc.6287, rm.m3$gsc.6157), sep="")), size=5, parse=TRUE)+
  theme(legend.position="none")

require(gridExtra)
pdf("~/MP_RStudio_May2014/Publication figures/GSC_reprod_BiotinElution.pdf", paper="a4r", width=20, height=10)
grid.arrange(plot1, plot2, plot3, ncol=3)
dev.off();

# Here are the corresponding correlations +++++++++++++++++++++++++++++++++++++++++
cor(m3$gsc.5832, m3$gsc.6157)
cor(m3$gsc.5832, m3$gsc.6287)
cor(m3$gsc.6157, m3$gsc.6287)

# Producing pdf ==========================================================================
pdf("~/MP_RStudio_May2014/CellSurfaceExtractionKit/BiotinElutionTryouts_gsc.pdf")
par(mfrow=c(1,3))
plot(log10(shared.5832$global.spectra.count))
plot(log10(shared.6157$global.spectra.count))
plot(log10(shared.6287$global.spectra.count))
dev.off();

