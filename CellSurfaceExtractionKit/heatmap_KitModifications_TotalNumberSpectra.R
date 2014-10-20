# Set your wd
setwd('/slow/home/khoermann/MP/CellSurfaceExtractionKit/'); getwd();

# Read in the data
P5824 <- read.csv2('/slow/home/khoermann/MP/CellSurfaceExtractionKit/P5824/20140325_Isobar_P5824.csv', header=T, stringsAsFactors=F)
ac.1 <- gsub("[,-].*",'',P5824$protein.acs)
gsc.1 <- P5824$global.spectra.count

P5826 <- read.csv2('/slow/home/khoermann/MP/CellSurfaceExtractionKit//P5826/20140325_Isobar_P5826.csv',header=T, stringsAsFactors=F)
ac.2 <- gsub("[,-].*",'',P5826$protein.acs)
gsc.2 <- P5826$global.spectra.count

P5828 <- read.csv2('/slow/home/khoermann/MP/CellSurfaceExtractionKit//P5828/20140325_Isobar_P5828.csv',header=T, stringsAsFactors=F)
ac.3 <- gsub("[,-].*",'',P5828$protein.acs)
gsc.3 <- P5828$global.spectra.count

P5830 <- read.csv2('/slow/home/khoermann/MP/CellSurfaceExtractionKit//P5830/20140325_Isobar_P5830.csv',header=T, stringsAsFactors=F)
ac.4 <- gsub("[,-].*",'',P5830$protein.acs)
gsc.4 <- P5830$global.spectra.count

P5832 <- read.csv2('/slow/home/khoermann/MP/CellSurfaceExtractionKit//P5832/20140325_Isobar_P5832.csv',header=T, stringsAsFactors=F)
ac.5 <- gsub("[,-].*",'',P5832$protein.acs)
gsc.5 <- P5832$global.spectra.count


gscs <- paste("gsc", 1:5, sep = ".")
gsc_all <- mget(gscs, envir = globalenv())  # Create a list containing gsc counts for every run
gsc_all <- lapply(gsc_all, as.numeric)

gsc_sum <- lapply(gsc_all, sum) # Calculate spectra # for each run
names(gsc_sum) <- c("P5824", "P5826", "P5828", "P5830", "P5832")
write.csv2(gsc_sum, file="GSC_sums_KitModifications.csv")

# This part is with too many GO Annotations giving no nice enrichment for the relatively few GO Terms
# associated with PM
# Loads GO Annotation===========================================================
sp2go.c <- read.csv2("/slow/home/khoermann/MP/Human_gocc.txt",check.names=F,sep="\t",stringsAsFactors=F)
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
pds <- c(1:5)
cns <- c("P5824", "P5826", "P5828", "P5830", "P5832")
go.pvals <- matrix(rep(1,5*length(go.descr)),ncol=5,dimnames=list(names(go.descr),cns))
acs <- paste("ac", 1:5, sep=".")
asc_all <- mget(acs, envir = globalenv())

# Fill the matrix
for (pd in pds){
    prot.set <- asc_all[[pd]]
    go <- goOneSet(prot.set)
    for (i in 1:dim(go)[1]) {
      go.pvals[go[i,"term"],pd] <- go[i,"p.val"]}
  }


# Save the matrix into a new variable and write it to a .csv file
hm.input <- go.pvals
colnames(hm.input) <- c("P5824", "P5826", "P5828", "P5830", "P5832")
write.csv2(hm.input, file="heatmap.GOCCs.KitModifications.csv")

# Continue working with the .csv file to save memory
hm.read <- read.csv2('heatmap.GOCCs.KitModifications.csv', header=TRUE, stringsAsFactors=F)

hm <- as.matrix(hm.read[,2:6, drop=FALSE])  # Take just the numerical columns

# Define the names for the columns extracted one row above
rownames(hm) <- hm.read[[1]]
colnames(hm) <- c("P5824", "P5826", "P5828", "P5830", "P5832")
best.pval <- apply(hm,1,min)  # Get the lowest p-value for each GO Term to see which cutoff makes sense
cut.off <- 0.0001
good <- best.pval<cut.off  # Check how many GO Terms pass the cutoff
sum(good)

useful <- hm[good,]  # Extract only the GO Terms that pass the filter
labels <- paste(rownames(useful),go.descr[rownames(useful)])  # Define the labels
heatmap(-log10(useful),scale="none",Colv=NA,col=rev(heat.colors(20)),labRow=labels)   # Create the heatmap

# Save the heatmap to a .pdf
pdf("heatmap_gocc_protocols.pdf",width=16,height=8,pointsize=8,useDingbats=F,)
heatmap(-log10(useful),scale="none",Colv=NA,col=rev(heat.colors(20)),labRow=labels)
dev.off()
