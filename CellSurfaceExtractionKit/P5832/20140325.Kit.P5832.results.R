# Process Samples obtained with Pierce Kit according to Florian's scheme
-------------------------------------------------------------------------------------

library(isobar)

#list.files of interest /MP/CellSurfaceExtractionKit/id.csv
# P5832: biotinylated + Biotin elution
setwd("/slow/home/khoermann/MP/CellSurfaceExtractionKit/id.csv/"); getwd(); 
id.files <- Sys.glob(c("P5832-[12].id.csv"))
id.names <- gsub(".*/(P[0-9]{4})-[12].*","\\1",id.files)
id.names

protein.group <- readProteinGroup(id.files)

require(RPostgreSQL)
proteinInfo(protein.group) <- getProteinInfoFromBioDb(protein.group,PostgreSQL(),user="biodbprod",password="#4biodbprod",
                                                      host="head1",dbname="cemmprod")
protein.groups <- list()
for (id.n in unique(id.names)) {
  protein.groups[[id.n]] <- readProteinGroup(id.files[id.names==id.n],template=protein.group)
}
protein.groups

protein.gs <- reporterProteins(protein.group)

result.df <- data.frame(protein.g=protein.gs,protein.acs=isobar:::.protein.acc(protein.gs,protein.group))
result.df <- cbind(result.df,proteinNameAndDescription(protein.group,protein.gs))
result.df$global.peptide.count <- peptide.count(protein.group,protein.gs,specificity=c(REPORTERSPECIFIC,GROUPSPECIFIC))
result.df$global.spectra.count <- spectra.count(protein.group,protein.gs,specificity=c(REPORTERSPECIFIC,GROUPSPECIFIC))
result.df$global.dNSAF <- calculate.dNSAF(protein.group)

ff <- list(peptide.count=function(...) peptide.count(...,specificity=c(REPORTERSPECIFIC,GROUPSPECIFIC)),
          spectra.count=function(...) spectra.count(...,specificity=c(REPORTERSPECIFIC,GROUPSPECIFIC)),
          dNSAF=function(...) calculate.dNSAF(...))

for (fff in names(ff)) {
  message(fff)
  my.f <- ff[[fff]]
  for (id.n in names(protein.groups)) {
    message(paste(fff,id.n,sep="."))
    result.df[[paste(fff,id.n,sep=".")]] <- my.f(protein.groups[[id.n]])[protein.gs]
  }
}

setwd("/slow/home/khoermann/MP/CellSurfaceExtractionKit/"); getwd(); 
write.csv2(result.df, file='20140325_Isobar_P5832.csv')

result.df[["protein.acs"]] <- gsub("[,-].*",'',result.df[['protein.acs']])
x <- result.df[['protein.acs']]

library(biomaRt)
unimart <- useMart("unimart",dataset="uniprot",
                   host="www.ebi.ac.uk",path="/uniprot/biomart/martservice")


# Take a look at the available attributes and filters
attr = listAttributes(unimart)
head(attr)
filters = listFilters(unimart)
head (filters)

go.terms=getBM(attributes=c("accession", "name", "go_name"),filters="accession",
                      values=x,mart=unimart)


# Organize the output as dataframe
DF <- as.data.frame(go.terms)

# Improve the representation to show one protein per row
DF2 <- aggregate(DF[3], DF[-3], 
                 FUN = function(X) paste(unique(X), collapse=", "))

write.csv2(DF2, file="20140325.Kit.P5832.GOterms.csv")

================================================================================
  
# Check how many of the identified proteins carried the Biotinylation modification
  
# Read in the .id.csv files
id.1 <- read.csv2("~/MP_RStudio_May2014/CellSurfaceExtractionKit/id.csv//P5832-1.id.csv", header=T, sep="\t", stringsAsFactors=F)
id.2 <- read.csv2("~/MP_RStudio_May2014/CellSurfaceExtractionKit/id.csv//P5832-2.id.csv", header=T, sep="\t", stringsAsFactors=F)

# Grep the biotinylated peptides
biot.1 <- id.1[grep("Biotin", id.1$modif), "accession"]
biot.2 <- id.2[grep("Biotin", id.2$modif), "accession"]
biot.1 <- gsub("[,-].*",'',biot.1) 
biot.2 <- gsub("[,-].*",'',biot.2)

total.biot <- sort(unique(c(biot.1, biot.2)))   # 81 proteins biotinylated (18% of PM proteins)
write.csv2(total.biot, "~/MP_RStudio_May2014/CellSurfaceExtractionKit/P5832/Biotinylated_ACs_P5832.csv")

