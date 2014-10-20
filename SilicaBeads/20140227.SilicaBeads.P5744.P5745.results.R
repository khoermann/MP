# Process Samples obtained with Silica Beads protocol according to Florian's scheme
-------------------------------------------------------------------------------------

library(isobar)

#list.files of interest /MP/Aminooxybiotin

getwd(); setwd("/slow/home/khoermann/MP/SilicaBeads/");
id.files <- Sys.glob(c("/slow/home/khoermann/MP/SilicaBeads/P5744-[12].id.csv","/slow/home/khoermann/MP/SilicaBeads/P5745-[12].id.csv"))
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

write.csv2(result.df, file='Isobar_P5744.P5745.csv')

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

go.terms.split <- ldply(strsplit(go.terms$go_name,","),function(x) {
  c(go.bp=paste0(sub("^P:","",x[grep("^P:",x)]),collapse=";"),
    go.cc=paste0(sub("^C:","",x[grep("^C:",x)]),collapse=";"),
    go.f=paste0(sub("^F:","",x[grep("^F:",x)]),collapse=";"))
})

go.cc=go.terms[grep("C:", go.terms$go_name),]

write.csv2(go.cc, file="20140227.SilicaBeads.P5744.P5745.GOCC.results.csv")

