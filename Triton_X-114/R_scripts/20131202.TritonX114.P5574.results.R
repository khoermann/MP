# Process Samples obtained with Leo's protocol (Triton X-114 phase separation) according to Florian's scheme
-------------------------------------------------------------------------------------

library(isobar)

#list.files of interest /MP/Triton_X-114

getwd(); setwd("/slow/home/khoermann/MP/Triton_X-114/");
id.files <- Sys.glob("/slow/home/khoermann/MP/Triton_X-114/P5574-[12].id.csv")
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

library(biomaRt)
unimart <- useMart("unimart",dataset="uniprot",
                   host="www.ebi.ac.uk",path="/uniprot/biomart/martservice")

result.df$protein.g <- as.character(result.df$protein.g)
go.terms <- ldply(strsplit(result.df$protein.g,","),function(x) {
  x <- unique(gsub("-.*","",x))
  message(paste0(x,collapse=","))
  c(go_name=paste0(getBM(attributes=c("go_name"),filters="accession",
                         values=x,mart=unimart),collapse=","),
    plasmid_name=paste0(getBM(attributes=c("plasmid_name"),filters="accession",
                               values=x,mart=unimart),collapse=","),
    organelle_name=paste0(getBM(attributes=c("organelle_name"),filters="accession",
                               values=x,mart=unimart),collapse=","))
})

go.terms.split <- ldply(strsplit(go.terms$go_name,","),function(x) {
  c(go.bp=paste0(sub("^P:","",x[grep("^P:",x)]),collapse=";"),
    go.cc=paste0(sub("^C:","",x[grep("^C:",x)]),collapse=";"),
    go.f=paste0(sub("^F:","",x[grep("^F:",x)]),collapse=";"))
})

result.df$unspecific.spectra.count <- spectra.count(protein.group,protein.gs,specificity=c(UNSPECIFIC))

head(go.terms.split)
writeIBSpectra(cbind(result.df,go.terms.split),file="20131202.TritonX114.P5574.results.csv",na="0")
system("/fast/bioinfo/sw/isobar/inst/pl/tab2xlsx.pl 20131202.TritonX114.P5574.results.csv")



