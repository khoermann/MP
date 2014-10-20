# Script to create a heatmap based on normalized spectral counts for the three GO categories: 
# PM, other membrane, other

# Read in the prepared data and extract the norm.sc column for every subset
# P5824
PM_P5824 <- read.csv2("/slow/home/khoermann/MP/CellSurfaceExtractionKit/P5824/PM_5824.csv", header=T, stringsAsFactors=F)
pm.1 <- PM_P5824$norm.sc
Non_PM_mem_P5824 <- read.csv2("/slow/home/khoermann/MP/CellSurfaceExtractionKit/P5824/Non_PM_mem_P5824.csv", header=T, stringsAsFactors=F)
n.1 <- Non_PM_mem_P5824$norm.sc
Other_P5824 <- read.csv2("/slow/home/khoermann/MP/CellSurfaceExtractionKit/P5824/Other_P5824.csv", header=T, stringsAsFactors=F)
o.1 <- Other_P5824$norm.sc

# P5826
PM_P5826 <- read.csv2("/slow/home/khoermann/MP/CellSurfaceExtractionKit/P5826//PM_5826.csv", header=T, stringsAsFactors=F)
pm.2 <- PM_P5826$norm.sc
Non_PM_mem_P5826 <- read.csv2("/slow/home/khoermann/MP/CellSurfaceExtractionKit/P5826//Non_PM_mem_P5826.csv", header=T, stringsAsFactors=F)
n.2 <- Non_PM_mem_P5826$norm.sc
Other_P5826 <- read.csv2("/slow/home/khoermann/MP/CellSurfaceExtractionKit/P5826//Other_P5826.csv", header=T, stringsAsFactors=F)
o.2 <- Other_P5826$norm.sc

# P5828
PM_P5828 <- read.csv2("/slow/home/khoermann/MP/CellSurfaceExtractionKit/P5828/PM_5828.csv", header=T, stringsAsFactors=F)
pm.3 <- PM_P5828$norm.sc
Non_PM_mem_P5828 <- read.csv2("/slow/home/khoermann/MP/CellSurfaceExtractionKit/P5828/Non_PM_mem_P5828.csv", header=T, stringsAsFactors=F)
n.3 <- Non_PM_mem_P5828$norm.sc
Other_P5828 <- read.csv2("/slow/home/khoermann/MP/CellSurfaceExtractionKit/P5828/Other_P5828.csv", header=T, stringsAsFactors=F)
o.3 <- Other_P5828$norm.sc

# P5830
PM_P5830 <- read.csv2("/slow/home/khoermann/MP/CellSurfaceExtractionKit/P5830/PM_5830.csv", header=T, stringsAsFactors=F)
pm.4 <- PM_P5830$norm.sc
Non_PM_mem_P5830 <- read.csv2("/slow/home/khoermann/MP/CellSurfaceExtractionKit/P5830/Non_PM_mem_P5830.csv", header=T, stringsAsFactors=F)
n.4 <- Non_PM_mem_P5830$norm.sc
Other_P5830 <- read.csv2("/slow/home/khoermann/MP/CellSurfaceExtractionKit/P5830/Other_P5830.csv", header=T, stringsAsFactors=F)
o.4 <- Other_P5830$norm.sc

# P5832
PM_P5832 <- read.csv2("/slow/home/khoermann/MP/CellSurfaceExtractionKit/P5832/PM_5832.csv", header=T, stringsAsFactors=F)
pm.5 <- PM_P5832$norm.sc
Non_PM_mem_P5832 <- read.csv2("/slow/home/khoermann/MP/CellSurfaceExtractionKit/P5832/Non_PM_mem_P5832.csv", header=T, stringsAsFactors=F)
n.5 <- Non_PM_mem_P5832$norm.sc
Other_P5832 <- read.csv2("/slow/home/khoermann/MP/CellSurfaceExtractionKit/P5832/Other_P5832.csv", header=T, stringsAsFactors=F)
o.5 <- Other_P5832$norm.sc


# Create lists containing the subsets according to GO annotation
pms <- paste("pm", 1:5, sep = ".")
pm_all <- mget(pms, envir = globalenv())

ns <- paste("n", 1:5, sep = ".")
n_all <- mget(ns, envir = globalenv())

os <- paste("o", 1:5, sep = ".")
o_all <- mget(os, envir = globalenv())

# Calculate the sum for all
pm_means <- lapply(pm_all, sum)
n_means <- lapply(n_all, sum)
o_means <- lapply(o_all, sum)

pds <- c("P5824", "P5826", "P5828", "P5830", "P5832")
exps <- c(1:5)
terms <- c("Plasma Membrane", "Non-PM membrane", "Others")
m.to.fill <- matrix(rep(0,5*length(terms)),ncol=5,dimnames=list(terms,pds))

for (exp in exps){ 
  m.to.fill[1, exp] <- pm_means[[exp]]
  m.to.fill[2, exp] <- n_means[[exp]]
  m.to.fill[3, exp] <- o_means[[exp]]
}

write.csv2(m.to.fill, file="/slow/home/khoermann/MP/CellSurfaceExtractionKit/matrixHeatmap_SumNormSC_KitModifications.csv")

png("hm_SumNormSC_KitModifications.png")
heatmap(m.to.fill,scale="none",Colv=NA,
        col=rev(heat.colors(70, alpha=.7)),
        labRow=terms, cexRow=1, cexCol=1, margins=c(5,10))   # Still needs some refinement in what tranformation to use on the data before creating the heatmap
dev.off();