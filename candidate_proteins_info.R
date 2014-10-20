require( RPostgreSQL )

candidate.genes <- read.table( '~/Add.Haplogen.July2013.txt', header = TRUE )

idb.conn <- dbConnect('PostgreSQL', host="head1", dbname="cemmprod", user="astukalov", password='#4astukalov' )

sql.query <- paste0( "
SELECT de.*
FROM biodbprod_20130121.dbentries AS de
WHERE de.primaryac IN (", paste( "'", candidate.genes$protein_ac, "'", sep='', collapse = "," ), ")" )

sql.query <- paste0( "
SELECT de.chromosome, COUNT( de.primaryac )
FROM biodbprod_20130121.dbentries AS de
GROUP BY de.chromosome" )

prot_info <- dbGetQuery( idb.conn, sql.query )

require( biomaRt )

ensembl.mart <- useMart( "ensembl", dataset =  'hsapiens_gene_ensembl' )

candate_genes.attribs <- getBM(
  attributes = c( 'uniprot_swissprot_accession', 'uniprot_genename', 'chromosome_name', 'description', 'entrezgene' ),
  filters = c( 'uniprot_swissprot_accession' ),
  values = candidate.genes$protein_ac, mart = ensembl.mart )

table( candate_genes.attribs$chromosome_name )

subset( candate_genes.attribs, chromosome_name == '15' )

listFilters( ensembl.mart )
bm.attrs <- listAttributes( ensembl.mart, what = c('name','description','page') )

grep( 'gene', listAttributes( ensembl.mart, what = c('name','description','page') )$name, value=TRUE )

str(listAttributes( ensembl.mart, showGroups=TRUE))

