## Analysis of data from Cell Surface extraction kit 

getwd(); setwd('/slow/home/khoermann/MP/')

# Read processed data
library(xlsx)

data <- read.xlsx2('20130705.mp.results.xlsx', sheetIndex=1, header=TRUE)

# Define IDs and ACCs
Kit.ids <- unique(as.character(data$ID))
Kit.acc <- as.character(data$Protein.Acs)

# Get annotation info
# Define the library
library(biomaRt)

# Define the dataset (Uniprot) and server to use
unimart <- useMart("unimart",dataset="uniprot",
                   host="www.ebi.ac.uk",path="/uniprot/biomart/martservice")

# Take a look at the available attributes and filters
attr = listAttributes(unimart)
head(attr)
filters = listFilters(unimart)
head (filters)

# getBM
results.Kit <- getBM(attributes=c('accession', 'name', 'gene_name', 'db2go_c__dm_description'), filters = 'accession', 
                 values = Kit.acc, mart=unimart)
as.data.frame(results.Kit)

mem.Kit <- results.Kit[grep('membrane', results.Kit$db2go_c__dm_description), ]
unique.mem.acc <- unique(mem.Kit$accession)

# Search results with a vector of GOCC annotations of interest
toMatch <- c('plasma membrane', 'cell surface', 'extracellular')
matches.Kit <- results.Kit[grep(paste(toMatch, collapse='|'), results.Kit$db2go_c__dm_description), ]
ids.matches.Kit <- sort(unique(matches.Kit$name))


# Filter the excel output file created by Florian's scheme for proteins with GOCC 'integral to membrane' and no other subcellular localiz.
itm.xlsx <- read.xlsx2('20131203_ITM_CellSurfaceKit.xlsx', sheetIndex=1, header=TRUE)
ids.itm.xlsx <- as.character(itm.xlsx$ID)
setdiff(ids.itm.xlsx, ids.matches.Kit)



--------------------------------------------------------
--------------------------------------------------------
  # Alternative way to extract all interesting GOCCs one by one
mem <- results[grep('membrane', results$db2go_c__dm_description), ]
unique.mem.ids <- unique(mem$name)

pm <- results[grep('plasma membrane', results$db2go_c__dm_description), ]
unique.pm.ids <- unique(pm$name)

cs <- results[grep('cell surface', results$db2go_c__dm_description), ]
unique.cs.ids <- unique(cs$name)

ext <- results[grep('extracellular', results$db2go_c__dm_description), ]
unique.ext.ids <- unique(ext$name)


all <- unique(c(unique.cs.ids, unique.ext.ids, unique.pm.ids))

-------------------------------------------------------
  
  

