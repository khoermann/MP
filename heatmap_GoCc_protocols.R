# Set your wd
setwd('/slow/home/khoermann/MP/'); getwd();

# Read in the data
require(xlsx)
P5134_6 <- read.xlsx2('/slow/home/khoermann/MP/CellSurfaceExtractionKit//P5134-5136/20130705.mp.results.xlsx', sheetIndex=1, stringsAsFactors=F)
P5134_6.acc <- gsub("[,-].*",'',P5134_6$Protein.Acs)
gsc.1 <- P5134_6$Global.Spectra.Count

P5826 <- read.csv2('/slow/home/khoermann/MP/CellSurfaceExtractionKit//P5826/20140325_Isobar_P5826.csv',header=T, stringsAsFactors=F)
P5826.acc <- gsub("[,-].*",'',P5826$protein.acs)
gsc.2 <- P5826$global.spectra.count

P5396_5455 <- read.xlsx2("/slow/home/khoermann/MP/Aminooxybiotin/P5396_P5455/20131202.Aminooxy.P5396.P5455.results.xlsx",sheetIndex=1, stringsAsFactors=F)
P5396_5455.acc <- gsub("[,-].*",'',P5396_5455$Protein.Acs)
gsc.3 <- P5396_5455$Global.Spectra.Count

P5767_8 <- read.csv2("/slow/home/khoermann/MP/Aminooxybiotin/P5767_P5768//201400324_Isobar_P5767.P5768.csv",header=T, stringsAsFactors=F)
P5767_8.acc <- gsub("[,-].*",'',P5767_8$protein.acs)
gsc.4 <- P5767_8$global.spectra.count

P5571 <- read.xlsx2("/slow/home/khoermann/MP/SilicaBeads//20131203.SilicaBeads.P5571.results.xlsx", sheetIndex=1, stringsAsFactors=F)
P5571.acc <- gsub("[,-].*",'',P5571$Protein.Acs)
gsc.5 <- P5571$Global.Spectra.Count

P5572 <- read.xlsx2("/slow/home/khoermann/MP/SilicaBeads//20131203.SilicaBeads.P5572.results.xlsx", sheetIndex=1, stringsAsFactors=F)
P5572.acc <- gsub("[,-].*",'',P5572$Protein.Acs)
gsc.6 <- P5572$Global.Spectra.Count

P5744 <- read.csv2("/slow/home/khoermann/MP/SilicaBeads/20140325_Isobar_P5744.csv", header=T, stringsAsFactors=F)
P5744.acc <- gsub("[,-].*",'',P5744$protein.acs)
gsc.7 <- P5744$global.spectra.count

P5745 <- read.csv2("/slow/home/khoermann/MP/SilicaBeads/20140325_Isobar_P5745.csv", header=T, stringsAsFactors=F)
P5745.acc <- gsub("[,-].*",'',P5745$protein.acs)
gsc.8 <- P5745$global.spectra.count
prot.index <- list(P5134_6.acc, P5826.acc, P5396_5455.acc, P5767_8.acc, P5571.acc, P5744.acc, P5572.acc, P5745.acc)

gscs <- paste("gsc", 1:8, sep = ".")
gsc_all <- mget(gscs, envir = globalenv())  # Create a list containing gsc counts for every run
gsc_all <- lapply(gsc_all, as.numeric)

gsc_sum <- lapply(gsc_all, sum) # Calculate spectra # for each run
names(gsc_sum) <- c("P5134_6", "P5826", "P5396_5455", "P5767_5768", "P5571", "P5572", "P5744", "P5745")
write.csv2(gsc_sum, file="GSC_sums.csv")


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


# Prepare a matrix 
pds <- c(1:8)
go.pvals <- matrix(rep(1,8*length(go.descr)),ncol=8,dimnames=list(names(go.descr),pds))

# Fill the matrix
for (pd in pds){
    prot.set <- prot.index[[pd]]
    go <- goOneSet(prot.set)
    for (i in 1:dim(go)[1]) {
      go.pvals[go[i,"term"],pd] <- go[i,"p.val"]}
  }


# Save the matrix into a new variable and write it to a .csv file
hm.input <- go.pvals
write.csv2(hm.input, file="heatmap.GOCCs.protocols.csv")

# Continue working with the .csv file to save memory
hm.read <- read.csv2('heatmap.GOCCs.protocols.csv', header=TRUE, stringsAsFactors=F)

hm <- as.matrix(hm.read[,2:9, drop=FALSE])  # Take just the numerical columns

# Define the names for the columns extracted one row above
rownames(hm) <- hm.read[[1]]
colnames(hm) <- c("P5134_6", "P5826", "P5396_5455", "P5767_8", "P5571", "P5744", "P5572", "P5745")
best.pval <- apply(hm,1,min)  # Get the lowest p-value for each GO Term to see which cutoff makes sense
cut.off <- 0.00001
good <- best.pval<cut.off  # Check how many GO Terms pass the cutoff
sum(good)

useful <- hm[good,]  # Extract only the GO Terms that pass the filter
labels <- paste(rownames(useful),go.descr[rownames(useful)])  # Define the labels
heatmap(-log10(useful),scale="none",Colv=NA,col=rev(heat.colors(20)),labRow=labels)   # Create the heatmap

# Save the heatmap to a .pdf
pdf("heatmap_gocc_protocols.pdf",width=16,height=8,pointsize=8,useDingbats=F,)
heatmap(-log10(useful),scale="none",Colv=NA,col=rev(heat.colors(20)),labRow=labels)
dev.off()
