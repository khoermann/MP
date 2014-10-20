# Looking at how many of the biotinylated proteins are PM
# Sort of control for whether the reagent is really non-membrane permeable as claimed by Thermo. 

# "classic" approach, 1st tryout, P5134
biotin.P5134 <- read.csv2("~/MP_RStudio_May2014/CellSurfaceExtractionKit/P5134-5136/P5134_BiotinylatedProteins.csv", header=T, stringsAsFactors=F)
PM.P5134 <- read.csv2("~/MP_RStudio_May2014/CellSurfaceExtractionKit/P5134-5136/PM_P5134.csv", header=T, stringsAsFactors=F)
int.P5134 <- intersect(biotin.P5134$x, PM.P5134$protein.acs)
diff.P5134 <- setdiff(biotin.P5134$x, PM.P5134$protein.acs)

# "classic" approach, 2nd tryout, P5826
biotin.P5826 <- read.csv2("~/MP_RStudio_May2014/CellSurfaceExtractionKit/P5826/Biotinylated_ACs_P5826.csv", header=T, stringsAsFactors=F)
PM.P5826 <- read.csv2("~/MP_RStudio_May2014/CellSurfaceExtractionKit/P5826/PM_5826.csv", header=T, stringsAsFactors=F)
int.P5826 <- intersect(biotin.P5826$x, PM.P5826$protein.acs)
diff.P5826 <- setdiff(biotin.P5826$x, PM.P5826$protein.acs)

# 1st tryout, P5832
biotin.P5832 <- read.csv2("~/MP_RStudio_May2014/CellSurfaceExtractionKit/P5832/Biotinylated_ACs_P5832.csv", header=T, stringsAsFactors=F)
PM.P5832 <- read.csv2("~/MP_RStudio_May2014/CellSurfaceExtractionKit/P5832/PM_5832.csv", header=T, stringsAsFactors=F)
int.P5832 <- intersect(biotin.P5832$x, PM.P5832$protein.acs)
diff.P5832 <- setdiff(biotin.P5832$x, PM.P5832$protein.acs)

# 2nd tryout, P6157
biotin.P6157 <- read.csv2("~/MP_RStudio_May2014/CellSurfaceExtractionKit/P6157/Biotinylated_ACs_P6157.csv", header=T, stringsAsFactors=F)
PM.P6157 <- read.csv2("~/MP_RStudio_May2014/CellSurfaceExtractionKit/P6157/PM_6157.csv", header=T, stringsAsFactors=F)
int.P6157 <- intersect(biotin.P6157$x, PM.P6157$protein.acs)
diff.P6157 <- setdiff(biotin.P6157$x, PM.P6157$protein.acs)

# 3rd tryout, P6287
biotin.P6287 <- read.csv2("~/MP_RStudio_May2014/CellSurfaceExtractionKit/P6287/Biotinylated_ACs_P6287.csv", header=T, stringsAsFactors=F)
PM.P6287 <- read.csv2("~/MP_RStudio_May2014/CellSurfaceExtractionKit/P6287/PM_P6287.csv", header=T, stringsAsFactors=F)
int.P6287 <- intersect(biotin.P6287$x, PM.P6287$protein.acs)
diff.P6287 <- setdiff(biotin.P6287$x, PM.P6287$protein.acs)

# Urea wash, P6284
biotin.P6284 <- read.csv2("~/MP_RStudio_May2014/CellSurfaceExtractionKit/P6284/Biotinylated_ACs_P6284.csv", header=T, stringsAsFactors=F)
PM.P6284 <- read.csv2("~/MP_RStudio_May2014/CellSurfaceExtractionKit/P6284/PM_P6284.csv", header=T, stringsAsFactors=F)
int.P6284 <- intersect(biotin.P6284$x, PM.P6284$protein.acs)
diff.P6284 <- setdiff(biotin.P6284$x, PM.P6284$protein.acs)

# Checking whether we see an enrichment for other cellular compartments
# Loads GO Annotation===========================================================
sp2go.c <- read.csv2("Human_gocc.txt",check.names=F,sep="\t",stringsAsFactors=F)
uc <- unique(sp2go.c[,3:4])
go.descr.c <- uc[[2]]
names(go.descr.c) <- uc[[1]]

head(sp2go.c)
sp2go <- sp2go.c[,c(1,3,4)]
head(sp2go)
go.descr <- go.descr.c


# Compares ======================================
# Calculation of p-values for each GO Term
goOneSet <- function(prots){
  go <- NULL
  for (term in names(go.descr)){
    inpw <- sp2go[sp2go[[2]]==term,1]
    pval <- phyper(q=length(intersect(prots,inpw))-1,m=length(inpw),n=20265-length(inpw),k=length(prots),lower.tail=F)
    go <- rbind(go,data.frame(term=term,descr=uc[grep(term, uc[,1]), 2], p.val=pval,num.inpw=length(inpw),
                              num.prots=length(prots),num.intersect=length(intersect(inpw,prots)),stringsAsFactors=F))
  }
  go
}

pdnr <- c("P5134", "P5826", "P5832", "P6157", "P6287")
pdnracs <- paste("diff", pdnr, sep=".")
lspdn <- mget(pdnracs, envir = globalenv())

gocc.pval <- c()
gocc.num <- c()

for (i in 1:length(pdnr)){
  df.acs <- lspdn[[i]]
  go <- goOneSet(df.acs)
  ord.go1 <- go[order(go$num.intersect, decreasing=T), ]
  ord.go2 <- go[order(go$p.val), ]
  top.ord.go1 <- ord.go1[1:6, ]
  top.ord.go2 <- ord.go2[1:6, ]
  top.ord.go1 -> gocc.num[[i]]
  top.ord.go2 -> gocc.pval[[i]]
}

write.table(gocc.num, file="~/MP_RStudio_May2014/CellSurfaceExtractionKit/BiotinylationModif_by_GoCC_num.txt")
write.table(gocc.pval, file="~/MP_RStudio_May2014/CellSurfaceExtractionKit/BiotinylationModif_by_GoCC_pval.txt")
