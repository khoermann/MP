## Getting information for KO candidates KBM7/Haplogen

# Read in the .csv file
data <- read.csv2("/slow/home/khoermann/MP/GenesSelected.Oct2013.csv", header=TRUE, sep=";")
data

# Create subsets
Haplogen <- subset(data, Letter.code=='H')
Haplogen.July <- subset(data, Letter.code=='H July')
Rabs <- subset(data, Letter.code=='R')

## Grep all row numbers where letter.code starts with a D
rows.D <- grep ('^D', data$Letter.code)
# Define the Diseases subset
Diseases <- data [rows.D,]

-----------------------------------------------------------
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

# Get localizations
results <- getBM(attributes=c('name', 'gene_name', 'db2go_c__dm_description'), filters = 'accession', values = data[,1], mart=unimart)
results

write.csv(results, file='20131028.localization.KBM7-Haplogen.cand.csv')

-----------------------------------------------------------------------------
  
# Change to Ensembl homo sapiens data base to get chromosome localization
ensembl.mart <- useMart( "ensembl", dataset =  'hsapiens_gene_ensembl' )

chr.results <- getBM(attributes = c('uniprot_genename', 'uniprot_swissprot_accession', 'chromosome_name'),
                     filters = c( 'uniprot_swissprot_accession' ),
                     values = data[,1], mart = ensembl.mart )

chr.results

write.csv(chr.results, file='20131028.chromosome_name.KBM7-Haplogen.cand.csv')  




