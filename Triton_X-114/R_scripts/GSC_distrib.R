# Script to compare GSC distribution for TX-114 protocol

setwd("~/MP_RStudio_May2014/Triton_X-114/"); getwd();

# Read in the data for the subgroups =====================================
# Start creating dataframes for each subgroup containing localiz and phase as factors
# PM
det.pm <- read.csv2("~/MP_RStudio_May2014/Triton_X-114//P5834_5835/PM_P5834.P5835.csv", header=T, stringsAsFactors=F)
aqu.pm <- read.csv2("~/MP_RStudio_May2014/Triton_X-114//P5836_5837/PM_P5836.P5837.csv", header=T, stringsAsFactors=F)
ctrl.pm <- read.csv2("~/MP_RStudio_May2014/Triton_X-114//P5838_5839/PM_P5838.P5839.csv", header=T, stringsAsFactors=F)

boxplot(log10(det.pm$global.spectra.count), 
        log10(aqu.pm$global.spectra.count), 
        log10(ctrl.pm$global.spectra.count), 
        ylab="Log10(Global Spectra Counts)", 
        names=c("det", "aqu", "ctrl"),
        col=c("firebrick1", "royalblue", "grey"), 
        main="Global Spectra Count distribution of PM fraction")

pm.df <- rbind(det.pm[, c("protein.acs", "global.spectra.count")], 
                aqu.pm[, c("protein.acs", "global.spectra.count")], 
                ctrl.pm[, c("protein.acs", "global.spectra.count")])
pm.df$localiz <- rep("PM", 2302)
pm.df$phase <- c(rep("det", 775), rep("aqu", 712), rep("ctrl", 815))

# PM only
det.pmo <- read.csv2("~/MP_RStudio_May2014/Triton_X-114/P5834_5835/PMOnly_P5834.P5835.csv", header=T, stringsAsFactors=F)
aqu.pmo <- read.csv2("~/MP_RStudio_May2014/Triton_X-114//P5836_5837/PMOnly_P5836.P5837.csv", header=T, stringsAsFactors=F)
ctrl.pmo <- read.csv2("~/MP_RStudio_May2014/Triton_X-114//P5838_5839/PMOnly_P5838.P5839.csv", header=T, stringsAsFactors=F)

pmo.df <- rbind(det.pmo[, c("protein.acs", "global.spectra.count")], 
               aqu.pmo[, c("protein.acs", "global.spectra.count")], 
               ctrl.pmo[, c("protein.acs", "global.spectra.count")])
pmo.df$localiz <- rep("PM Only", 1092)
pmo.df$phase <- c(rep("det", 474), rep("aqu", 269), rep("ctrl", 349))


# Membrane
det.mem <- read.csv2("~/MP_RStudio_May2014/Triton_X-114//P5834_5835/Mem_P5834.P5835.csv", header=T, stringsAsFactors=F)
aqu.mem <- read.csv2("~/MP_RStudio_May2014/Triton_X-114//P5836_5837/Mem_P5836.P5837.csv", header=T, stringsAsFactors=F)
ctrl.mem <- read.csv2("~/MP_RStudio_May2014/Triton_X-114//P5838_5839/Mem_P5838.P5839.csv", header=T, stringsAsFactors=F)

boxplot(log10(det.mem$global.spectra.count), 
        log10(aqu.mem$global.spectra.count), 
        log10(ctrl.mem$global.spectra.count), 
        ylab="Log10(Global Spectra Counts)", 
        names=c("det", "aqu", "ctrl"),
        col=c("firebrick1", "royalblue", "grey"), 
        main="Global Spectra Count distribution of Membrane fraction")

mem.df <- rbind(det.mem[, c("protein.acs", "global.spectra.count")], 
               aqu.mem[, c("protein.acs", "global.spectra.count")], 
               ctrl.mem[, c("protein.acs", "global.spectra.count")])
mem.df$localiz <- rep("Membrane", 2631)
mem.df$phase <- c(rep("det", 1245), rep("aqu", 585), rep("ctrl", 801))

# Non-PM-membrane
det.npm <- read.csv2("~/MP_RStudio_May2014/Triton_X-114//P5834_5835/Non_PM_mem_P5834.P5835.csv", header=T, stringsAsFactors=F)
aqu.npm <- read.csv2("~/MP_RStudio_May2014/Triton_X-114//P5836_5837/Non_PM_mem_P5836.P5837.csv", header=T, stringsAsFactors=F)
ctrl.npm <- read.csv2("~/MP_RStudio_May2014/Triton_X-114//P5838_5839//Non_PM_mem_P5838_5839.csv", header=T, stringsAsFactors=F)

boxplot(log10(det.npm$global.spectra.count), 
        log10(aqu.npm$global.spectra.count), 
        log10(ctrl.npm$global.spectra.count), 
        ylab="Log10(Global Spectra Counts)", 
        names=c("det", "aqu", "ctrl"),
        col=c("firebrick1", "royalblue", "grey"), 
        main="Global Spectra Count distribution of non-PM membrane fraction")

npm.df <- rbind(det.npm[, c("protein.acs", "global.spectra.count")], 
               aqu.npm[, c("protein.acs", "global.spectra.count")], 
               ctrl.npm[, c("protein.acs", "global.spectra.count")])
npm.df$localiz <- rep("Membrane-PM", 1268)
npm.df$phase <- c(rep("det", 675), rep("aqu", 238), rep("ctrl", 355))

# Other
det.other <- read.csv2("~/MP_RStudio_May2014/Triton_X-114//P5834_5835/Other_P5834.P5835.csv", header=T, stringsAsFactors=F)
aqu.other <- read.csv2("~/MP_RStudio_May2014/Triton_X-114//P5836_5837//Other_P5836.P5837.csv", header=T, stringsAsFactors=F)
ctrl.other <- read.csv2("~/MP_RStudio_May2014/Triton_X-114//P5838_5839//Other_P5838_5839.csv", header=T, stringsAsFactors=F)

boxplot(log10(det.other$global.spectra.count), 
        log10(aqu.other$global.spectra.count), 
        log10(ctrl.other$global.spectra.count), 
        ylab="Log10(Global Spectra Counts)", 
        names=c("det", "aqu", "ctrl"),
        col=c("firebrick1", "royalblue", "grey"), 
        main="Global Spectra Count distribution of non-membrane proteins")

other.df <- rbind(det.other[, c("protein.acs", "global.spectra.count")], 
               aqu.other[, c("protein.acs", "global.spectra.count")], 
               ctrl.other[, c("protein.acs", "global.spectra.count")])
other.df$localiz <- rep("Non-membrane", 2081)
other.df$phase <- c(rep("det", 431), rep("aqu", 836), rep("ctrl", 814))

# SLCs/ABCs
det.int <- read.csv2("~/MP_RStudio_May2014/Triton_X-114//P5834_5835/SLCsABCs_P5834.P5835.csv", header=T, stringsAsFactors=F)
aqu.int <- read.csv2("~/MP_RStudio_May2014/Triton_X-114//P5836_5837//SLCsABCs_P5836.P5837.csv", header=T, stringsAsFactors=F)
ctrl.int <- read.csv2("~/MP_RStudio_May2014/Triton_X-114//P5838_5839//SLCsABCs_P5838.P5839.csv", header=T, stringsAsFactors=F)

boxplot(log10(det.int$global.spectra.count), 
        log10(aqu.int$global.spectra.count), 
        log10(ctrl.int$global.spectra.count), 
        ylab="Log10(Global Spectra Counts)", 
        names=c("det", "aqu", "ctrl"),
        col=c("firebrick1", "royalblue", "grey"), 
        main="Global Spectra Count distribution of SLC/ABC proteins")

int.df <- rbind(det.int[, c("protein.acs", "global.spectra.count")], 
               aqu.int[, c("protein.acs", "global.spectra.count")], 
               ctrl.int[, c("protein.acs", "global.spectra.count")])
int.df$localiz <- rep("SLC/ABC", 113)
int.df$phase <- c(rep("det", 82), rep("aqu", 11), rep("ctrl", 20))      

# Bind to one big dataframe
all.df <- rbind(pm.df, mem.df, npm.df, other.df, int.df, pmo.df)

# Plot via ggplot =========================================================
require("MASS")
library("ggplot2")

# By localization
plot <- ggplot(all.df, aes(x = localiz, y = log10(global.spectra.count))) +
  geom_boxplot(aes(fill = phase)) +
  theme_bw() +
  xlab("") +
  theme(axis.text=element_text(size=8)) +
  scale_fill_manual(values = c("blue", "grey", "red"))

pdf("~/MP_RStudio_May2014/Triton_X-114/pdf/GSC_allLocalizPhases_inclPMOnly.pdf")
plot
dev.off();

# By phase
plot2 <- ggplot(all.df, aes(x = phase, y = log10(global.spectra.count))) +
  geom_boxplot(aes(fill = localiz)) +
  theme_bw() +
  xlab("")
  
pdf("~/MP_RStudio_May2014/Triton_X-114/pdf/GSC_localizByPhases.pdf")
plot2
dev.off();
  
 
# Look at GSC distrib per phase =========================================
require("RColorBrewer")
# Detergent phase
pdf("~/MP_RStudio_May2014/Triton_X-114/pdf/GSC_distrib_det.pdf", width=9, height=9)
boxplot(log10(det.mem$global.spectra.count), 
        log10(det.pm$global.spectra.count), 
        log10(det.npm$global.spectra.count), 
        log10(det.other$global.spectra.count), 
        log10(det.int$global.spectra.count),
        ylab="Log10(Global Spectra Counts)", 
        names=c("Membrane", "PM", "Membrane-PM", "Non-membrane", "SLC/ABC"),
        col=brewer.pal(5, "Reds"), 
        main="Global Spectra Count distribution of detergent phase")
dev.off();

# Aquous phase ==========================================================
pdf("~/MP_RStudio_May2014/Triton_X-114/pdf/GSC_distrib_aqu.pdf", width=9, height=9)
boxplot(log10(aqu.mem$global.spectra.count), 
        log10(aqu.pm$global.spectra.count), 
        log10(aqu.npm$global.spectra.count), 
        log10(aqu.other$global.spectra.count), 
        log10(aqu.int$global.spectra.count),
        ylab="Log10(Global Spectra Counts)", 
        ylim=c(0,3),
        names=c("Membrane", "PM", "Membrane-PM", "Non-membrane", "SLC/ABC"),
        col=brewer.pal(5, "Blues"), 
        main="Global Spectra Count distribution of aquous phase")
dev.off();

# Ctrl ==================================================================
pdf("~/MP_RStudio_May2014/Triton_X-114/pdf/GSC_distrib_ctrl.pdf", width=9, height=9)
boxplot(log10(ctrl.mem$global.spectra.count), 
        log10(ctrl.pm$global.spectra.count), 
        log10(ctrl.npm$global.spectra.count), 
        log10(ctrl.other$global.spectra.count), 
        log10(ctrl.int$global.spectra.count),
        ylab="Log10(Global Spectra Counts)", 
        ylim=c(0,3),
        names=c("Membrane", "PM", "Membrane-PM", "Non-membrane", "SLC/ABC"),
        col=brewer.pal(5, "Greys"), 
        main="Global Spectra Count distribution of control")
dev.off();