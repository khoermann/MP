# Identify 'poster childs' = PM protein identified with purification methods, but not with standard procedure

# Read in all the data
P5134_5136 <- read.csv2("/slow/home/khoermann/MP/CellSurfaceExtractionKit/P5134-5136/PM_P5134_5136.csv", header=T, stringsAsFactors=F)
P5826 <- read.csv2("/slow/home/khoermann/MP/CellSurfaceExtractionKit/P5826/PM_5826.csv", header=T, stringsAsFactors=F)
P4569 <- read.csv2("/slow/home/khoermann/MP/PM_P4569_kbm7.csv", header=T, stringsAsFactors=F)

# Identify the "poster childs" and the proteins that might be particularly intersting
length(setdiff(P5134_5136$Protein.Acs, P4569$AC))
new.PM <- P5134_5136[grep(paste(setdiff(P5134_5136$Protein.Acs, P4569$AC), collapse="|"), P5134_5136$Protein.Acs), ]
new.slc <- new.PM[grep("SLC", new.PM$Gene), ]
new.abc <- new.PM[grep("ABC", new.PM$Gene), ]
new.kinases <- new.PM[grep("kinase", new.PM$Description), ]
write.csv2(rbind(new.slc, new.abc, new.kinases), file="/slow/home/khoermann/MP/PCsOfInterest_1stTryout.csv")

length(setdiff(P5826$protein.acs, P4569$AC))
add <- P5826[grep(paste(setdiff(P5826$protein.acs, P4569$AC), collapse="|"), P5826$protein.acs), ]



