## Analysis of tryout with Silica bead protocol

getwd(); setwd('/slow/home/khoermann/MP/SilicaBeads/')

# Read processed data
library(xlsx)

data <- read.xlsx2('20131203.SilicaBeads.P5570.results.xlsx', sheetIndex=1, header=TRUE)

# Define IDs and ACCs
P5570.ids <- unique(as.character(data$ID))
P5570.acc <- as.character(data$Protein.Acs)

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
results.P5570 <- getBM(attributes=c('name', 'gene_name', 'db2go_c__dm_description'), filters = 'accession', 
                 values = P5570.acc, mart=unimart)
as.data.frame(results.P5570)

mem.P5570 <- results.P5570[grep('membrane', results.P5570$db2go_c__dm_description), ]
unique.mem.ids <- unique(mem.P5570$name)

# Search results with a vector of GOCC annotations of interest
toMatch <- c('plasma membrane', 'cell surface', 'extracellular')
matches.P5570 <- results.P5570[grep(paste(toMatch, collapse='|'), results.P5570$db2go_c__dm_description), ]
ids.matches.P5570 <- sort(unique(matches.P5570$name))

# Extend to the term 'integral to membrane' and then manually check back whether there is no further subcellular localiz. known
itm.P5570 <- results.P5570[grep('integral to membrane', results.P5570$db2go_c__dm_description), ]

# Limit down to unique IDs
itm.P5570.ids <- unique(itm.P5570$name)

# Identify those not already covered by GOCC annotations of interest
cand <- setdiff(itm.P5570.ids, ids.matches.P5570)

# Filter the excel output file created by Florian's scheme 
itm.xlsx <- read.xlsx2('20131203_ITM_P5570X114.xlsx', sheetIndex=1, header=TRUE)
ids.itm.xlsx <- as.character(itm.xlsx$ID)
setdiff(ids.itm.xlsx, ids.matches.P5570)



--------------------------------------------------------
--------------------------------------------------------
  # Alternative way to extract all interesting GOCCs one by one
mem.amino <- results[grep('membrane', results$db2go_c__dm_description), ]
unique.mem.ids <- unique(mem.amino$gene_name)

pm.amino <- results[grep('plasma membrane', results$db2go_c__dm_description), ]
unique.pm.ids <- unique(pm.amino$gene_name)

cs.amino <- results[grep('cell surface', results$db2go_c__dm_description), ]
unique.cs.ids <- unique(cs.amino$gene_name)

ext.amino <- results[grep('extracellular', results$db2go_c__dm_description), ]
unique.ext.ids <- unique(ext.amino$gene_name)

int.amino <- results[grep('integral to membrane', results$db2go_c__dm_description), ]
unique.int.ids <- unique(int.amino$gene_name)

all <- unique(c(unique.cs.ids, unique.ext.ids, unique.int.ids, unique.pm.ids))

-------------------------------------------------------
  
  

