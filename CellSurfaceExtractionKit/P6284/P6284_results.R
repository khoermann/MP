# Process Samples obtained with Pierce Kit according to Florian's scheme
-------------------------------------------------------------------------------------

library(isobar)

#list.files of interest /MP/CellSurfaceExtractionKit/id.csv
# P6284: 40 Mio wt cells labeled, Cell Surface Extraction Kit, washing with 8M Urea, Biotin elution
setwd("~/MP_RStudio_May2014/CellSurfaceExtractionKit//id.csv"); getwd(); 
id.files <- Sys.glob(c("P6284-[12].id.csv"))
id.names <- gsub(".*/(P[0-9]{4})-[12].*","\\1",id.files)
id.names

protein.group <- readProteinGroup(id.files)

require(RPostgreSQL)
proteinInfo(protein.group) <- getProteinInfoFromBioDb(protein.group,PostgreSQL(),user="biodbprod",password="#4biodbprod",
                                                      host="msprod",dbname="cemmprod")
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

dim(result.df)  # 646 proteins detected
setwd("~/MP_RStudio_May2014/CellSurfaceExtractionKit/P6284/"); getwd(); 
write.csv2(result.df, file='Isobar_P6284.csv')

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

go.terms=getBM(attributes=c("accession", "db2go_c__dm_description"),filters="accession",
               values=x,mart=unimart)

# Organize the output as dataframe
DF <- as.data.frame(go.terms)

# Improve the representation to show one protein per row
DF2 <- aggregate(DF[2], DF[-2], 
                 FUN = function(X) paste(unique(X), collapse=", "))   
colnames(DF2) <- c("protein.acs", "Go.Cc")

combi <- merge(result.df, DF2, by="protein.acs")     # Combine Isobar data and annotations into one dataframe
dim(combi)
write.csv2(combi, file="P6284.GOterms.csv")       # Contains 637 entries

--------------------------------------------------------------------------
  
# Check how many of the identified proteins carried the Biotinylation modification
  
# Read in the .id.csv files
id.1 <- read.csv2("~/MP_RStudio_May2014/CellSurfaceExtractionKit/id.csv//P6284-1.id.csv", header=T, sep="\t", stringsAsFactors=F)
id.2 <- read.csv2("~/MP_RStudio_May2014/CellSurfaceExtractionKit/id.csv//P6284-2.id.csv", header=T, sep="\t", stringsAsFactors=F)

# Grep the biotinylated peptides
biot.1 <- id.1[grep("Biotin", id.1$modif), "accession"]
biot.2 <- id.2[grep("Biotin", id.2$modif), "accession"]
biot.1 <- gsub("[,-].*",'',biot.1) 
biot.2 <- gsub("[,-].*",'',biot.2)

total.biot <- sort(unique(c(biot.1, biot.2)))   # 38 proteins biotinylated (10% of PM proteins)
write.csv2(total.biot, "~/MP_RStudio_May2014/CellSurfaceExtractionKit/P6284/Biotinylated_ACs_P6284.csv")
