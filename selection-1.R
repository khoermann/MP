library(SparseM)
library(RBGL)
library(Rgraphviz)

source("../net-r/expandGraph.R")

base.name.net <- "/home/colinge/databases/interact-public-human-withcomplexes-2012-07-24"
load(file=paste(base.name.net,"-graph.data",sep=""))
n.alldb.G <- nodes(alldb.G)
load(paste(base.name.net,"-transmat.data",sep=""))

# ID mapping ===============================

bcrabl <- "A9UF07"
sp2name <- read.csv("/home/colinge/databases/sp2name-biomart-sponly-2012-11-26.txt",check.names=F,sep="\t",stringsAsFactors=F)
sp.names <- sp2name[["Gene name"]]
names(sp.names) <- sp2name[["Accession"]]
sp.names[sp.names==""] <- sp2name[["Accession"]][sp.names==""]
sp.names[bcrabl] <- "BCR-ABL"
names.sp <- sp2name[["Accession"]]
names(names.sp) <- sp2name[["Gene name"]]
names.sp["BCR-ABL"] <- bcrabl
sp.names["Q9Y3C1"] <- "NOP16"
id.sp <- sp2name[["Accession"]]
names(id.sp) <- sp2name[["Entry name"]]

# number of KEGG pathways in which a protein is found =================
sp2kegg <- read.csv("/home/colinge/databases/uniprot_to_kegg_human-2012-08-23.csv",check.names=F,sep="\t",stringsAsFactors=F)
sp.kegg <- sp2kegg[!is.na(sp2kegg[[2]]),]
for (i in 1:dim(sp.kegg)[1]){
  l <- nchar(sp.kegg[i,2])
  sp.kegg[i,2] <- paste(paste(rep("0",5-l),collapse=""),sp.kegg[i,2],sep="")
}
kegg2names <- read.csv("/home/colinge/databases/kegg_names-2012-08-23.txt",check.names=F,sep="\t",stringsAsFactors=F)
for (i in 1:dim(kegg2names)[1]){
  l <- nchar(kegg2names[i,1])
  kegg2names[i,1] <- paste(paste(rep("0",5-l),collapse=""),kegg2names[i,1],sep="")
}
kegg.names <- kegg2names[[2]]
names(kegg.names) <- kegg2names[[1]]

# Average number of pathways hit
kids <- unique(sp.kegg[[2]])
kpw <- NULL
for (k in kids)
  kpw <- c(kpw,list(sp.kegg[sp.kegg[[2]]==k,1]))
names(kpw) <- kids
num.in.pw <- table(unlist(kpw))
num.in.pw <- num.in.pw[num.in.pw>0]
summary(num.in.pw)

# Candidate genes =============================

cand <- read.csv("datasets/intersect.candidates.20130423.kbm7.txt",check.names=F,header=F,sep="\t",stringsAsFactors=F)[[1]]
length(cand)

# no nucleus-associated protein
nucleus <- union(read.csv("datasets/nucleus-proteins.txt",check.names=F,sep="\t",stringsAsFactors=F)[[2]],read.csv("datasets/nuclear-membrane.txt",check.names=F,sep="\t",stringsAsFactors=F)[[2]])
cand.non.nucl <- setdiff(cand,nucleus)
length(cand.non.nucl)
num.pw.cand <- num.in.pw[intersect(cand.non.nucl,names(num.in.pw))]
num.pw.cand
sp.names[names(num.pw.cand)[num.pw.cand>1]]

# Haplogen collection
haplo <- read.csv("datasets/picked-clones-2013-04-17-converted.txt",check.names=F,header=F,sep="\t",stringsAsFactors=F)[[1]]
cand.haplo <- intersect(cand.non.nucl,haplo)
length(cand.haplo)
cand.remain <- setdiff(cand.non.nucl,haplo)
length(cand.remain)

# Decoy complexes with 1 Haplogen clone at least
corum <- read.csv("/home/colinge/databases/corum-feb-2012.csv",check.names=F,sep="\t",stringsAsFactors=F)
decoy.cor <- NULL
decoy.cor.names <- NULL
for (i in 1:dim(corum)[1]){
  cplx <- unlist(strsplit(corum[i,5],","))
  if (length(intersect(cand.haplo,cplx)) > 0){
    decoy.cor <- c(decoy.cor,cplx)
    decoy.cor.names <- c(decoy.cor.names,corum[i,2])
  }
}
print(decoy.cor.names)
length(intersect(decoy.cor,cand.remain))
cand.remain <- setdiff(cand.remain,decoy.cor)
length(cand.remain)

# Decoy proteins in direct interaction with 1 Haplogen clone
deg <- degree(alldb.G)
max.deg <- quantile(deg,prob=0.95)
cand.haplo.G <- intersect(intersect(n.alldb.G,cand.haplo),names(deg[deg<max.deg]))
decoy.ppi <- unique(unlist(adj(alldb.G,cand.haplo.G)))
length(intersect(cand.remain,decoy.ppi))
cand.remain <- setdiff(cand.remain,decoy.ppi)
length(cand.remain)

# Add one family of interest avoiding redundancy in complexes and PPI
family <- read.csv("datasets/Rab.within.sorted.candidates.20130423",check.names=F,header=F,sep="\t",stringsAsFactors=F)[[1]]
family <- setdiff(family,nucleus)
length(family)
#family <- read.csv("datasets/Sec.within.sorted.candidates.20130423",check.names=F,header=F,sep="\t",stringsAsFactors=F)[[1]]
#length(family)
#intersect(family,haplo)
fam.haplo <- intersect(family,haplo)
fam.haplo

fam.cor <- NULL
fam.cor.names <- NULL
for (i in 1:dim(corum)[1]){
  cplx <- unlist(strsplit(corum[i,5],","))
  if (length(intersect(family,cplx)) > 0){
    fam.cor <- c(decoy.cor,cplx)
    fam.cor.names <- c(decoy.cor.names,corum[i,2])
  }
}
fam.red <- union(fam.haplo,intersect(family,fam.cor))
fam.red
fam.haplo.G <- intersect(intersect(n.alldb.G,fam.haplo),names(deg[deg<max.deg]))
fam.ppi <- unique(unlist(adj(alldb.G,fam.haplo.G)))
fam.red <- union(fam.haplo,intersect(family,fam.ppi))
fam.red

fam.to.add <- setdiff(family,fam.red)
fam.to.add
sel.fam <- intersect(fam.to.add,n.alldb.G)

fam.added <- sparseSelectionCorum(fam.to.add,5)
fam.added <- c(fam.added,sparseSelectionPPI(setdiff(fam.to.add,fam.added),15-length(fam.added)))
sel.fam <- union(fam.haplo,fam.added)
sel.fam

# Remaining
remaining <- setdiff(cand.remain,family)
remain.added <- sparseSelectionCorum(remaining,5)
remain.added <- c(remain.added,sparseSelectionPPI(setdiff(remaining,remain.added),5))
sel.remain <- remain.added
sel.remain

# Final selection
sel.cand <- unique(c(cand.haplo,sel.fam,sel.remain))
sel.cand
sp.names[sel.cand]

num.in.pw[intersect(sel.cand,names(num.in.pw))]
sp.names[names(num.in.pw[intersect(sel.cand,names(num.in.pw))])]



three <- names(num.pw.cand)[num.pw.cand>3]
common.pws <- NULL
for (i in 1:length(kpw)){
  if (sum(three %in% kpw[[i]])>0)
    common.pws <- union(common.pws,kegg.names[names(kpw)[i]])
}

three <- names(num.pw.cand)[num.pw.cand>=5]
common.pws <- NULL
for (i in 1:length(kpw)){
  if (sum(three %in% kpw[[i]])>0)
    common.pws <- union(common.pws,kegg.names[names(kpw)[i]])
}

# Greedy selection of candidates ------------------------------------------------
sparseSelectionCorum <- function(nodes,num=1){

  pos <- size <- NULL
  for (i in 1:dim(corum)[1]){
    cplx <- unlist(strsplit(corum[i,5],","))
    if (length(intersect(nodes,cplx)) > 0){
      pos <- c(pos,i)
      size <- c(size,length(intersect(cplx,nodes)))
    }
  }
  o <- order(size,decreasing=T)
  selection <- NULL
  i <- 1
  while ((length(selection) < num) && (i <= dim(corum)[1]){
    cplx <- unlist(strsplit(corum[o[i],5],","))
    selection <- union(selection,intersect(cplx,nodes)[1])
    i <- i+1
  }
  selection[!is.na(selection)]
}


sparseSelectionPPI <- function(nodes,num=1){

  nodes <- intersect(nodes,n.alldb.G)
  print(nodes)
  a <- adj(alldb.G,nodes)
  #print(a)
  pos <- size <- NULL
  for (i in 1:length(a)){
    pos <- c(pos,i)
    size <- c(size,length(intersect(a[[i]],nodes)))
  }
  print(size)
  o <- order(size,decreasing=T)
  print(o)
  selection <- NULL
  i <- 1
  while ((length(selection) < num) && (i <= length(size))){
    selection <- union(selection,intersect(a[[o[i]]],nodes)[1])
    i <- i+1
  }
  selection[!is.na(selection)]
}

