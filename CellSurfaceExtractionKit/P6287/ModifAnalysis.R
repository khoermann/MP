# Check how many of the identified proteins carried the Biotinylation modification

# Read in the .id.csv files
id.1 <- read.csv2("~/MP_RStudio_May2014/CellSurfaceExtractionKit/id.csv//P5134-1.id.csv", header=T, sep="\t", stringsAsFactors=F)
id.2 <- read.csv2("~/MP_RStudio_May2014/CellSurfaceExtractionKit/id.csv//P5134-2.id.csv", header=T, sep="\t", stringsAsFactors=F)

# Grep the biotinylated peptides
biot.1 <- id.1[grep("Biotin", id.1$modif), "accession"]
biot.2 <- id.2[grep("Biotin", id.2$modif), "accession"]
biot.1 <- gsub("[,-].*",'',biot.1) 
biot.2 <- gsub("[,-].*",'',biot.2)

total.biot <- sort(unique(c(biot.1, biot.2))) # Yields 85 biotinyl. proteins (22% of PM proteins)
write.csv2(total.biot, "~/MP_RStudio_May2014/CellSurfaceExtractionKit/P5134-5136/P5134_BiotinylatedProteins.csv") 
