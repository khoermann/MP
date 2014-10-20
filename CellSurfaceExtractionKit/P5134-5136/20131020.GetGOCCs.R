# Define Uniprot IDs
acc <- read.csv("/slow/home/khoermann/Acc.trafficking.Rabs.HAP1.csv", header=FALSE)

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
results <- getBM(attributes=c('name', 'gene_name', 'db2go_c__dm_description'), filters = 'accession', values = acc, mart=unimart)
results
as.data.frame(results)

write.csv(results, file='20131020.gocc.trafficking.Rab.HAP1.csv')






