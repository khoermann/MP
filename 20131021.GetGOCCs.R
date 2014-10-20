# Define Uniprot IDs
acc <- read.csv('/home/khoermann/Acc.to.get.csv', header=FALSE)
acc

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
results <- getBM(attributes=c('name', 'gene_name', 'protein_name'), filters = 'accession', values = acc, mart=unimart)
results

write.csv(results, file='20131021.protein_name.csv')

# Change to Ensembl homo sapiens data base to get chromosome localization
ensembl.mart <- useMart( "ensembl", dataset =  'hsapiens_gene_ensembl' )

chr.results <- getBM(attributes = c('uniprot_genename', 'uniprot_swissprot_accession', 'chromosome_name'),
  filters = c( 'uniprot_swissprot_accession' ),
  values = acc, mart = ensembl.mart )

chr.results

write.csv(chr.results, file='20131021.chromosome_name.csv')

# Take a look at the available attributes and filters for the ensembl database
attr.ensembl = listAttributes(ensembl.mart)
head(attr.ensembl)
filters.ensembl = listFilters(ensembl.mart)
head (filters.ensembl)







