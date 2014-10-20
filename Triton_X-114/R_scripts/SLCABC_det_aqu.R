setwd("/slow/home/khoermann/MP/Triton_X-114/"); getwd();

slc.det <- read.csv2("SLCsABCs_P5834.P5835.csv", header=T, stringsAsFactors=F)  # 82 proteins
slc.aqu <- read.csv2("SLCsABCs_P5836.P5837.csv", header=T, stringsAsFactors=F)  # 11 proteins
slc.10 <- read.csv2("SLCsABCs_P5574.csv", header=T, stringsAsFactors=F)   # 40 proteins
P4569 <- read.csv2("/slow/home/khoermann/MP/PM_P4569_kbm7.csv", header=T, stringsAsFactors=F)

length(intersect(slc.10$Protein.Acs, slc.det$protein.acs))  # All 40 from first tryout det phase are retrieved in 2nd tryout det phase
both <- slc.det[grep(paste(intersect(slc.det$protein.acs, slc.aqu$protein.acs), collapse="|"), slc.det$protein.acs), ]
dim(slc.det); dim(slc.aqu); dim(both); dim(slc.10)


