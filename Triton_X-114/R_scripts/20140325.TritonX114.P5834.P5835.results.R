# Process Samples obtained with Triton X-114 protocol according to Florian's scheme
-------------------------------------------------------------------------------------

library(isobar)

#list.files of interest /MP/Triton_X-114/id.csv
# P5834/P5835: detergent phases
setwd("/slow/home/khoermann/MP/Triton_X-114/id.csv"); getwd(); 
id.files <- Sys.glob(c("P5834-[12].id.csv", "P5835-[12].id.csv"))
id.names <- gsub(".*/(P[0-9]{4})-[34].*","\\1",id.files)
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

setwd("/slow/home/khoermann/MP/Triton_X-114/"); getwd(); 
write.csv2(result.df, file='20140325_Isobar_P5834.P5835.csv')


