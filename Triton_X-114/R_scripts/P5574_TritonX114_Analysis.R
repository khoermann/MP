## Analysis of tryout with Leo's protocol (Triton X-114 phase separation)
# 10 Mio cells, detergent phase

getwd(); setwd('/slow/home/khoermann/MP/Triton_X-114/')

# Read processed data
library(xlsx)

triton.df <- read.xlsx2('20131202.TritonX114.P5574.results.xlsx', sheetIndex=1, header=TRUE)
str(triton.df)
triton.df$Protein.Acs <- gsub("[,-].*",'',triton.df$Protein.Acs); dim(triton.df);

# Look for SLCs and ABCs
length(unique(grep("SLC", triton.df$Gene)))   # 36 SLCs
length(unique(grep("ABC", triton.df$Gene)))   # 4 ABCs
slc.5574 <- triton.df[grep("SLC", triton.df$Gene),]
abc.5574 <- triton.df[grep("ABC", triton.df$Gene),]
write.csv2(rbind(slc.5574, abc.5574), file="SLCsABCs_P5574.csv")  

mem <- triton.df[grep('membrane', triton.df$Go.Cc), ]; dim(mem);    # 733 proteins
write.csv2(mem, file="Mem_P5574.csv")
