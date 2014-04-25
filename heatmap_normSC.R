# Script to create a heatmap based on normalized spectral counts for the three GO categories: PM, other membrane, other

# Read in the prepared data and extract the norm.sc column for every subset
# P5134-6
PM_P5134_5136 <- read.csv2("/slow/home/khoermann/MP/CellSurfaceExtractionKit/P5134-5136/PM_P5134_5136.csv", header=T, stringsAsFactors=F)
pm.1 <- PM_P5134_5136$norm.sc
Non_PM_mem_P5134_5136 <- read.csv2("/slow/home/khoermann/MP/CellSurfaceExtractionKit/P5134-5136/Non_PM_mem_P5134_5136.csv", header=T, stringsAsFactors=F)
n.1 <- Non_PM_mem_P5134_5136$norm.sc
Other_P5134_5136 <- read.csv2("/slow/home/khoermann/MP/CellSurfaceExtractionKit/P5134-5136/Other_P5134_5136.csv", header=T, stringsAsFactors=F)
o.1 <- Other_P5134_5136$norm.sc

# P5826
PM_P5826 <- read.csv2("/slow/home/khoermann/MP/CellSurfaceExtractionKit/P5826//PM_5826.csv", header=T, stringsAsFactors=F)
pm.2 <- PM_P5826$norm.sc
Non_PM_mem_P5826 <- read.csv2("/slow/home/khoermann/MP/CellSurfaceExtractionKit/P5826//Non_PM_mem_P5826.csv", header=T, stringsAsFactors=F)
n.2 <- Non_PM_mem_P5826$norm.sc
Other_P5826 <- read.csv2("/slow/home/khoermann/MP/CellSurfaceExtractionKit/P5826//Other_P5826.csv", header=T, stringsAsFactors=F)
o.2 <- Other_P5826$norm.sc

# P5396_P5455
PM_P5396_5455 <- read.csv2("/slow/home/khoermann/MP/Aminooxybiotin/P5396_P5455/PM_5396.P5455.csv", header=T, stringsAsFactors=F)
pm.3 <- PM_P5396_5455$norm.sc
Non_PM_mem_P5396_5455 <- read.csv2("/slow/home/khoermann/MP/Aminooxybiotin/P5396_P5455/Non_PM_mem_P5396.P5455.csv", header=T, stringsAsFactors=F)
n.3 <- Non_PM_mem_P5396_5455$norm.sc
Other_P5396_5455 <- read.csv2("/slow/home/khoermann/MP/Aminooxybiotin/P5396_P5455/Other_P5396.P5455.csv", header=T, stringsAsFactors=F)
o.3 <- Other_P5396_5455$norm.sc

# P5767_P5768
PM_P5767_5768 <- read.csv2("/slow/home/khoermann/MP/Aminooxybiotin/P5767_P5768/PM_P5767.P5768.csv", header=T, stringsAsFactors=F)
pm.4 <- PM_P5767_5768$norm.sc
Non_PM_mem_P5767_5768 <- read.csv2("/slow/home/khoermann/MP/Aminooxybiotin/P5767_P5768/Non_PM_mem_P5767.P5768.csv", header=T, stringsAsFactors=F)
n.4 <- Non_PM_mem_P5767_5768$norm.sc
Other_P5767_5768 <- read.csv2("/slow/home/khoermann/MP/Aminooxybiotin/P5767_P5768/Other_P5767.P5768.csv", header=T, stringsAsFactors=F)
o.4 <- Other_P5767_5768$norm.sc

# P5571
PM_P5571 <- read.csv2("/slow/home/khoermann/MP/SilicaBeads/PM_P5571.csv", header=T, stringsAsFactors=F)
pm.5 <- PM_P5571$norm.sc
Non_PM_mem_P5571 <- read.csv2("/slow/home/khoermann/MP/SilicaBeads/Non_PM_mem_P5571.csv", header=T, stringsAsFactors=F)
n.5 <- Non_PM_mem_P5571$norm.sc
Other_P5571 <- read.csv2("/slow/home/khoermann/MP/SilicaBeads/Other_P5571.csv", header=T, stringsAsFactors=F)
o.5 <- Other_P5571$norm.sc

# P5572
PM_P5572 <- read.csv2("/slow/home/khoermann/MP/SilicaBeads/PM_P5572.csv", header=T, stringsAsFactors=F)
pm.6 <- PM_P5572$norm.sc
Non_PM_mem_P5572 <- read.csv2("/slow/home/khoermann/MP/SilicaBeads/Non_PM_mem_P5572.csv", header=T, stringsAsFactors=F)
n.6 <- Non_PM_mem_P5572$norm.sc
Other_P5572 <- read.csv2("/slow/home/khoermann/MP/SilicaBeads/Other_P5572.csv", header=T, stringsAsFactors=F)
o.6 <- Other_P5572$norm.sc

# P5744
PM_P5744 <- read.csv2("/slow/home/khoermann/MP/SilicaBeads/PM_P5744.csv", header=T, stringsAsFactors=F)
pm.7 <- PM_P5744$norm.sc
Non_PM_mem_P5744 <- read.csv2("/slow/home/khoermann/MP/SilicaBeads/Non_PM_mem_P5744.csv", header=T, stringsAsFactors=F)
n.7 <- Non_PM_mem_P5744$norm.sc
Other_P5744 <- read.csv2("/slow/home/khoermann/MP/SilicaBeads/Other_P5744.csv", header=T, stringsAsFactors=F)
o.7 <- Other_P5744$norm.sc

# P5745
PM_P5745 <- read.csv2("/slow/home/khoermann/MP/SilicaBeads/PM_P5745.csv", header=T, stringsAsFactors=F)
pm.8 <- PM_P5745$norm.sc 
Non_PM_mem_P5745 <- read.csv2("/slow/home/khoermann/MP/SilicaBeads/Non_PM_mem_P5745.csv", header=T, stringsAsFactors=F)
n.8 <- Non_PM_mem_P5745$norm.sc
Other_P5745 <- read.csv2("/slow/home/khoermann/MP/SilicaBeads/Other_P5745.csv", header=T, stringsAsFactors=F)
o.8 <- Other_P5745$norm.sc

# Create lists containing the subsets according to GO annotation
pms <- paste("pm", 1:8, sep = ".")
pm_all <- mget(pms, envir = globalenv())

ns <- paste("n", 1:8, sep = ".")
n_all <- mget(ns, envir = globalenv())

os <- paste("o", 1:8, sep = ".")
o_all <- mget(os, envir = globalenv())

# Calculate the sum for all
pm_means <- lapply(pm_all, sum)
n_means <- lapply(n_all, sum)
o_means <- lapply(o_all, sum)

pds <- c("P5134_5136", "P5826", "P5396_5455", "P5767_5768", "P5571", "P5572", "P5744", "P5745")
exps <- c(1:8)
terms <- c("Plasma Membrane", "Non-PM membrane", "Others")
m.to.fill <- matrix(rep(0,8*length(terms)),ncol=8,dimnames=list(terms,pds))

for (exp in exps){ 
  m.to.fill[1, exp] <- pm_means[[exp]]
  m.to.fill[2, exp] <- n_means[[exp]]
  m.to.fill[3, exp] <- o_means[[exp]]
}

write.csv2(m.to.fill, file="/slow/home/khoermann/MP/matrixHeatmap_SumNormSC.csv")

heatmap(m.to.fill,scale="none",Colv=NA,col=rev(heat.colors(70, alpha=.7)),labRow=terms, cexRow=1, cexCol=1)   # Still needs some refinement in what tranformation to use on the data before creating the heatmap
