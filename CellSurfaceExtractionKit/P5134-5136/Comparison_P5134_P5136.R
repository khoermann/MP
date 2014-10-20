## Script to compare proteins identified and their abundance between P5134 and P5136

library(xlsx)

# Read in the .xlsx files
P5134.xlsx <- read.xlsx2 ('~/MP_RStudio_May2014/CellSurfaceExtractionKit/P5134-5136/P5134-1-P5134-2-1.xlsx', header=TRUE, sheetName='SampleCompare', sheetIndex=1, stringsAsFactors=F)
P5136.xlsx <- read.xlsx2 ('~/MP_RStudio_May2014/CellSurfaceExtractionKit/P5134-5136//P5136-1-P5136-2.xlsx', header=TRUE, sheetName='SampleCompare', sheetIndex=1)

# Defining dataframes with ACCs
df.P5134 <- data.frame(P5134.xlsx$AC)
df.P5136 <- data.frame(P5136.xlsx$AC)

---------------------------------------------------
  ---------------------------------------------
# Find proteins that are unique to the two datasets
# Find proteins exclusive to P5134
rows.in.P5134.that.are.not.in.P5136  <- function(df.P5134,df.P5136)
{
  P5134.vec <- apply(df.P5134, 1, paste, collapse = "")
  P5136.vec <- apply(df.P5136, 1, paste, collapse = "")
  P5134.without.P5136.rows <- df.P5134[!P5134.vec %in% P5136.vec,]
  return(P5134.without.P5136.rows)
}
rows.in.P5134.that.are.not.in.P5136(df.P5134,df.P5136)

# Create a csv file with results
write.csv(rows.in.P5134.that.are.not.in.P5136(df.P5134,df.P5136), file='unique.P5134.ids.csv')

----------------------------------------
# Find proteins exclusive to P5136
rows.in.P5136.that.are.not.in.P5134  <- function(df.P5134,df.P5136)
{
  P5134.vec <- apply(df.P5134, 1, paste, collapse = "")
  P5136.vec <- apply(df.P5136, 1, paste, collapse = "")
  P5136.without.P5134.rows <- df.P5136[!P5136.vec %in% P5134.vec,]
  return(P5136.without.P5134.rows)
}
rows.in.P5136.that.are.not.in.P5134(df.P5134,df.P5136)

write.csv(rows.in.P5136.that.are.not.in.P5134(df.P5134,df.P5136), file='unique.P5136.ids.csv')

------------------------------------------------
# Find proteins that are common in the two datasets

common.rows  <- function(df.P5134,df.P5136)
{
  P5134.vec <- apply(df.P5134, 1, paste, collapse = "")
  P5136.vec <- apply(df.P5136, 1, paste, collapse = "")
  P5134.and.P5136.rows <- df.P5134[P5134.vec %in% P5136.vec,]
  return(P5134.and.P5136.rows)
}
common.rows(df.P5134,df.P5136)

write.csv(common.rows(df.P5134,df.P5136), file='/slow/home/khoermann/MP/common.acc.P5134.P5136.csv')

------------------------------------------------
  
# Get GOCCs
common <- read.csv2('/slow/home/khoermann/MP/common.acc.P5134.P5136.csv', sep=',', header=TRUE)
common.acc <- as.character(common[,2])

P5134.only <- read.csv2('/slow/home/khoermann/MP/unique.P5134.ids.csv', sep=',', header=TRUE)
P5134.only.acc <- as.character(P5134.only[,2])

P5136.only <- read.csv2('/slow/home/khoermann/MP/unique.P5136.ids.csv', sep=',', header=TRUE)
P5136.only.acc <- as.character(P5136.only[,2])

total <- c(common.acc, P5134.only.acc, P5136.only.acc)


# Get annotation info
# Define the library
library(biomaRt)

# Define the dataset (Uniprot) and server to use
unimart <- useMart("unimart",dataset="uniprot",
                   host="www.ebi.ac.uk",path="/uniprot/biomart/martservice")

results <- getBM(attributes=c('name', 'gene_name', 'db2go_c__dm_description'), filters = 'accession', 
                 values = total, mart=unimart)

search.pm <- results[grep ('plasma membrane', results$db2go_c__dm_description), ]
head(search.pm)

search.cs <- results[grep ('cell surface', results$db2go_c__dm_description), ]
head(search.cs)

search.ex <- results[grep ('extracellular', results$db2go_c__dm_description), ]
head(search.ex)

search.int <- results[grep ('integral to membrane', results$db2go_c__dm_description), ]
head(search.int)

unique(c(search.pm$name, search.cs$name, search.ex$name, search.int$name))
