# Checking the number of PM proteins in KBM7 proteome data
# Testing different ways/approaches to do that
setwd('~/MP_RStudio_May2014/'); getwd();

# Read in the KBM7 data
kbm7.data <- read.csv2('P4569_kbm7.csv', header=TRUE, stringsAsFactors=F)   # Dataset has 10351 protein entries (9870 unique Uniprot AC)
# Process and isolate unique ACs
kbm7.data$AC <- gsub("[,-].*",'',kbm7.data$AC)  

----------------------------
# 1st option: get GO Annotations for all proteins in the KBM7 proteome
library(biomaRt)
unimart <- useMart("unimart",dataset="uniprot",
                   host="www.ebi.ac.uk",path="/uniprot/biomart/martservice")
results.kbm7 <- getBM(attributes=c("accession", "db2go_c__dm_description"), filters = 'accession', 
                       values = kbm7.data$AC, mart=unimart)
pre.df <- as.data.frame(results.kbm7)

# Improve the representation to show one protein per row
kbm7.go.df <- aggregate(pre.df[2], pre.df[-2], 
                 FUN = function(X) paste(unique(X), collapse=", "))   
# Not GO.CC annotated proteins are missed! 8893 proteins with GO.CC annotations
colnames(kbm7.go.df) <- c("AC", "Go.Cc")

combi <- merge(kbm7.data, kbm7.go.df, by="AC")     # Combine Isobar data and annotations into one dataframe
write.csv2(combi, file="P4569_kbm7.GOterms.csv")

# Use the file created above =====================================================
annot.kbm7 <- read.csv2("P4569_kbm7.GOterms.csv", header=T, stringsAsFactors=F)

PM.only <- c('plasma membrane', 'cell surface', 'cell membrane')
matches.pm.only <- annot.kbm7[grep(paste(PM.only, collapse='|'), annot.kbm7$Go.Cc), ]   # 1676 matches
write.csv2(matches.pm.only, file="PMOnly_P4569_kbm7.csv")

mem <- annot.kbm7[grep('membrane', annot.kbm7$Go.Cc), ]   # 3734 proteins

# Search results with a vector of GOCC annotations of interest
toMatch <- c('plasma membrane', 'cell surface', 'extracellular', 'cell membrane')
matches.pm <- annot.kbm7[grep(paste(toMatch, collapse='|'), annot.kbm7$Go.Cc), ]   # 2729 matches
length(setdiff(mem$AC, matches.pm$AC)) # 1696 non-PM membrane proteins

annot.kbm7$matched_mask <- sapply( strsplit(annot.kbm7$Go.Cc, ';', fixed = TRUE ),
                              function ( Ccs ) length( Ccs ) > 0 && all( Ccs %in% 'integral to membrane' ) )
kbm7.itm <- subset(annot.kbm7, matched_mask)   # No matches for this subset

write.csv2(matches.pm, file="PM_P4569_kbm7.csv")

# Search results with a vector of GOCC annotations of interest
toMatch <- c('plasma membrane', 'cell surface', 'extracellular', 'cell membrane')
matches.pm <- kbm7.go.df[grep(paste(toMatch, collapse='|'), kbm7.go.df$db2go_c__dm_description), ]
length(unique(matches.pm[['accession']]))     # Yields 2592 proteins

# Limiting it to strictly terms that contain 'plasma membrane'
only.pm <- kbm7.go.df[grep('plasma membrane', kbm7.go.df$db2go_c__dm_description),]   # Leaves 1525

# Identify membrane proteins
matches.mem <- kbm7.go.df[grep('membrane', kbm7.go.df$db2go_c__dm_description),]   # 3552 proteins

# Calculate percentages for subsets above
message("Percentage of strictly PM:"); dim(only.pm)[1]/length(kbm7.acc)
message("Percentage of PM +extracellular +cell surface:");dim(matches.pm)[1]/length(kbm7.acc)
message("Percentage of membrane:"); dim(matches.mem)[1]/length(kbm7.acc)


------------------------------------
#########################################
# 2nd option: Go to Swissprot and query HUMAN + subcellular localization 'plasma membrane'  
#########################################
# Read in the downloaded lists
sp.pm <- read.table('Swissprot_PM_human_ACC_20140327.txt', stringsAsFactors=F)    # Contains 2845 Protein ACs (27.03.2014)
sp.pmcsext <- read.table('Swissprot_PMCSEXT_human_ACC_20140328.txt', stringsAsFactors=F)    # Contains 6418 Protein ACs (28.03.2014)

# Overlap with KBM7 proteome
message('In KBM7 proteome:'); length(intersect(kbm7.acc, sp.pm$V1)); # Gives 900 proteins
message ('Percentage:'); length(intersect(kbm7.acc, sp.pm$V1))/length(kbm7.acc)   # Corresponding to 9.1%

message('In KBM7 proteome:'); length(intersect(kbm7.acc, sp.pmcsext$V1)); # Gives ONLY 971 proteins, even though #of ACs is 2.26x larger
message ('Percentage:'); length(intersect(kbm7.acc, sp.pmcsext$V1))/length(kbm7.acc)   # Corresponding to 9.8%

# Check the overlap between 1st and 2nd option
length(intersect(only.pm$accession, sp.pm$V1))    # Only gives 786 proteins 

# Double check: subject the downloaded lists to biomaRt analysis
library(biomaRt)
unimart <- useMart("unimart",dataset="uniprot",
                   host="www.ebi.ac.uk",path="/uniprot/biomart/martservice")
results.pm <- getBM(attributes=c("accession", "name", "db2go_c__dm_description"), filters = 'accession', 
                      values = sp.pm$V1, mart=unimart)
df <- as.data.frame(results.pm)
check <- aggregate(df[3], df[-3], 
                        FUN = function(X) paste(unique(X), collapse=", "))  

# Identify those not having a PM annotation
not.pm <- check[grep('plasma membrane', check$db2go_c__dm_description, invert=T),]  # 286 
length(grep(paste(toMatch, collapse='|'), not.pm$db2go_c__dm_description))  # Connected to PM: 67

-----------------------------------------------
  
# Extract SLCs:
length(grep("SLC", kbm7.data$Gene_symbol))  # 144 SLCs (1.39%)
length(grep("ABC", kbm7.data$Gene_symbol))  # 27 ABCs (0.26%)
slc.subset <- kbm7.data[grep("SLC", kbm7.data$Gene_symbol), ]
write.csv2(slc.subset, file="SLCs.kbm7.csv")
