library(xlsx)

# Define dataset to get GOCCs for
# Start with those proteins being at least 2-fold more abundant in extracted samples
fold2_sc <- read.xlsx2("/slow/home/khoermann/MP/fold2_sc_P5375-1-P5375-2-P5134-1-P5134-2-P5136-1-P5136-2.only_selected.xlsx", 
                       header=TRUE, sheetIndex=1, sheetName='SampleCompare')
acc <- as.character(fold2_sc$AC)

# Full list
comp <- read.xlsx2("/slow/home/khoermann/MP/P5375-1-P5375-2-P5134-1-P5134-2-P5136-1-P5136-2-1.xlsx", 
header=TRUE, sheetIndex=1, sheetName='SampleCompare')
acc.comp <- as.character(comp$AC)

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
results <- getBM(attributes=c('name', 'gene_name', 'db2go_c__dm_description'), filters = 'accession', 
                 values = acc, mart=unimart)
as.data.frame(results) -> localiz.2fold_sc

index.pm <- which(localiz.2fold_sc$db2go_c__dm_description=='plasma membrane')
index.int.mem <- which(localiz.2fold_sc$db2go_c__dm_description=='integral to membrane')
index.mem <- which(localiz.2fold_sc$db2go_c__dm_description=='membrane')
index.cell.surface <- which(localiz.2fold_sc$db2go_c__dm_description=='cell surface')

unique(localiz.2fold_sc$gene_name[c(index.cell.surface, index.pm, index.int.mem)])


results.comp <- getBM(attributes=c('name', 'gene_name', 'db2go_c__dm_description'), filters = 'accession', 
                 values = acc.comp, mart=unimart)
as.data.frame(results.comp) -> localiz.comp

index.pm.comp <- which(localiz.comp$db2go_c__dm_description=='plasma membrane')
index.int.mem.comp <- which(localiz.comp$db2go_c__dm_description=='integral to membrane')
index.mem.comp <- which(localiz.comp$db2go_c__dm_description=='membrane')
index.cell.surface.comp <- which(localiz.comp$db2go_c__dm_description=='cell surface')

unique(localiz.comp$gene_name[c(index.cell.surface.comp, index.pm.comp, index.int.mem.comp)])







