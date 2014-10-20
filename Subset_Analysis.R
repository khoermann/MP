# Identify subsets within the extracted PM fraction

# First tryout of Pierce Kit, SDS elution =============================================================================
# Then for the GPI anchored, connect to Biomart and get the Go.Mf annotations. ! Annotations seem to be poor here! 
library(biomaRt)
unimart <- useMart("unimart",dataset="uniprot",
                   host="www.ebi.ac.uk",path="/uniprot/biomart/martservice")

# Take a look at the available attributes and filters
attr = listAttributes(unimart)
head(attr)
filters = listFilters(unimart)
head (filters)

go.terms=getBM(attributes=c("accession", "db2go_f__dm_description"),filters="accession",
               values=pm.5826$Protein.Acs,mart=unimart)

# Organize the output as dataframe
DF <- as.data.frame(go.terms)

# Improve the representation to show one protein per row
DF2 <- aggregate(DF[2], DF[-2], 
                 FUN = function(X) paste(unique(X), collapse=", "))   
colnames(DF2) <- c("Protein.Acs", "Go.Mf")
compl <- merge(pm.5134, DF2, by="Protein.Acs")
grep("GPI anchor", compl$Go.Mf)
write.csv2(compl, file="~/MP_RStudio_May2014/CellSurfaceExtractionKit/P5134-5136/PM_P5134_5136_Mf.csv")

# If you consider P5134 separately ======================================================================
pm.5134 <- read.csv2("~/MP_RStudio_May2014/CellSurfaceExtractionKit/P5134-5136/PM_P5134.csv", stringsAsFactors=F, header=T)
pmonly.5134 <- read.csv2("~/MP_RStudio_May2014/CellSurfaceExtractionKit/P5134-5136/PMOnly_P5134.csv", stringsAsFactors=F, header=T)

length(grep("integral", pmonly.5134$Go.Cc))
not.integral <- pmonly.5134[grep("integral", pmonly.5134$Go.Cc, invert=T), ]
to.assign <- setdiff(pm.5134$protein.acs, pmonly.5134$protein.acs)
new.df <- pm.5134[grep(paste(to.assign, collapse="|"), pm.5134$protein.acs), ]
length(grep("extracellular", new.df$Go.Cc))
length(grep("cell surface", new.df$Go.Cc))
length(grep("glycoprotein", new.df$Description))

go.terms=getBM(attributes=c("accession", "db2go_f__dm_description"),filters="accession",
               values=pm.5134$protein.acs,mart=unimart)

# Organize the output as dataframe
DF <- as.data.frame(go.terms)

# Improve the representation to show one protein per row
DF2 <- aggregate(DF[2], DF[-2], 
                 FUN = function(X) paste(unique(X), collapse=", "))   
colnames(DF2) <- c("protein.acs", "Go.Mf")
compl <- merge(pm.5134, DF2, by="protein.acs")
grep("GPI anchor", compl$Go.Mf)
write.csv2(compl, file="~/MP_RStudio_May2014/CellSurfaceExtractionKit/P5134-5136/PM_P5134_Mf.csv")


# Second tryout of Pierce Kit, SDS elution =============================================================================
# P5826
pm.5826 <- read.csv2("~/MP_RStudio_May2014/CellSurfaceExtractionKit/P5826/PM_5826.csv", stringsAsFactors=F, header=T)
pmonly.5826 <- read.csv2("~/MP_RStudio_May2014/CellSurfaceExtractionKit/P5826/PMOnly_P5826.csv", stringsAsFactors=F, header=T)

length(grep("integral", pmonly.5826$Go.Cc))
to.assign <- setdiff(pm.5826$protein.acs, pmonly.5826$protein.acs)
new.df <- pm.5826[grep(paste(to.assign, collapse="|"), pm.5826$protein.acs), ]
length(grep("extracellular", new.df$Go.Cc))
length(grep("cell surface", new.df$Go.Cc))
length(grep("glycoprotein", new.df$Description))
go.terms=getBM(attributes=c("accession", "db2go_f__dm_description"),filters="accession",
               values=pm.5826$protein.acs,mart=unimart)

# Organize the output as dataframe
DF <- as.data.frame(go.terms)

# Improve the representation to show one protein per row
DF2 <- aggregate(DF[2], DF[-2], 
                 FUN = function(X) paste(unique(X), collapse=", "))   
colnames(DF2) <- c("protein.acs", "Go.Mf")
compl <- merge(pm.5826, DF2, by="protein.acs")
grep("GPI anchor", compl$Go.Mf)
write.csv2(compl, file="~/MP_RStudio_May2014/CellSurfaceExtractionKit/P5826/PM_P5826_Mf.csv")


# First tryout of Aminooxybiotin ======================================================================================
# P5396/P5455
pm.amino2 <- read.csv2("~/MP_RStudio_May2014/Aminooxybiotin/P5396_P5455/PM_5396.P5455.csv", stringsAsFactors=F, header=T)

length(grep("integral", pm.amino1$Go.Cc))
length(grep("extracellular", pm.amino1$Go.Cc))
length(grep("cell surface", pm.amino1$Go.Cc))
length(grep("glycoprotein", pm.amino1$Description))

go.terms=getBM(attributes=c("accession", "db2go_f__dm_description"),filters="accession",
               values=pm.amino1$Protein.Acs,mart=unimart)

# Organize the output as dataframe
DF <- as.data.frame(go.terms)

# Improve the representation to show one protein per row
DF2 <- aggregate(DF[2], DF[-2], 
                 FUN = function(X) paste(unique(X), collapse=", "))   
colnames(DF2) <- c("Protein.Acs", "Go.Mf")
compl <- merge(pm.amino1, DF2, by="Protein.Acs")
grep("GPI", compl$Go.Mf)
write.csv2(compl, file="~/MP_RStudio_May2014/Aminooxybiotin/P5396_P5455/PM_5396_P5455_Mf.csv")

# P5455 alone ==========================================================================
pm.5455 <- read.csv2("~/MP_RStudio_May2014/Aminooxybiotin/P5396_P5455/P5455/PM_P5455.csv", stringsAsFactors=F, header=T)
pmonly.5455 <- read.csv2("~/MP_RStudio_May2014/Aminooxybiotin/P5396_P5455/P5455/PMOnly_P5455.csv", stringsAsFactors=F, header=T)

length(grep("integral", pmonly.5455$Go.Cc))
to.assign <- setdiff(pm.5455$protein.acs, pmonly.5455$protein.acs)
new.df <- pm.5455[grep(paste(to.assign, collapse="|"), pm.5455$protein.acs), ]
length(grep("extracellular", new.df$Go.Cc))
length(grep("cell surface", new.df$Go.Cc))
length(grep("glycoprotein", new.df$Description))

go.terms=getBM(attributes=c("accession", "db2go_f__dm_description"),filters="accession",
               values=pm.5455$protein.acs,mart=unimart)

# Organize the output as dataframe
DF <- as.data.frame(go.terms)

# Improve the representation to show one protein per row
DF2 <- aggregate(DF[2], DF[-2], 
                 FUN = function(X) paste(unique(X), collapse=", "))   
colnames(DF2) <- c("protein.acs", "Go.Mf")
compl <- merge(pm.5455, DF2, by="protein.acs")
grep("GPI", compl$Go.Mf)
write.csv2(compl, file="~/MP_RStudio_May2014/Aminooxybiotin/P5396_P5455/P5455/PM_P5455_Mf.csv")


# Second tryout of Aminooxybiotin
# P5767_P5768
pm.amino2 <- read.csv2("~/MP_RStudio_May2014/Aminooxybiotin/P5767_P5768/PM_P5767.P5768-34.csv", stringsAsFactors=F, header=T)

length(grep("integral", pm.amino2$Go.Cc))
length(grep("extracellular", pm.amino2$Go.Cc))
length(grep("cell surface", pm.amino2$Go.Cc))
length(grep("glycoprotein", pm.amino2$Description))

go.terms=getBM(attributes=c("accession", "db2go_f__dm_description"),filters="accession",
               values=pm.amino2$protein.acs,mart=unimart)

# Organize the output as dataframe
DF <- as.data.frame(go.terms)

# Improve the representation to show one protein per row
DF2 <- aggregate(DF[2], DF[-2], 
                 FUN = function(X) paste(unique(X), collapse=", "))   
colnames(DF2) <- c("protein.acs", "Go.Mf")
compl <- merge(pm.amino2, DF2, by="protein.acs")
grep("GPI", compl$Go.Mf)
write.csv2(compl, file="~/MP_RStudio_May2014/Aminooxybiotin/P5767_P5768/PM_P5767.P5768-34_Mf.csv")

# P5767 alone =============================================================
pm.5767 <- read.csv2("~/MP_RStudio_May2014/Aminooxybiotin/P5767_P5768/P5767/PM_P5767.csv", stringsAsFactors=F, header=T)
pmonly.5767 <- read.csv2("~/MP_RStudio_May2014/Aminooxybiotin/P5767_P5768/P5767/PMOnly_P5767.csv", stringsAsFactors=F, header=T)

length(grep("integral", pmonly.5767$Go.Cc))
to.assign <- setdiff(pm.5767$protein.acs, pmonly.5767$protein.acs)
new.df <- pm.5767[grep(paste(to.assign, collapse="|"), pm.5767$protein.acs), ]
length(grep("extracellular", new.df$Go.Cc))
length(grep("cell surface", new.df$Go.Cc))
length(grep("glycoprotein", new.df$Description))

go.terms=getBM(attributes=c("accession", "db2go_f__dm_description"),filters="accession",
               values=pm.5767$protein.acs,mart=unimart)

# Organize the output as dataframe
DF <- as.data.frame(go.terms)

# Improve the representation to show one protein per row
DF2 <- aggregate(DF[2], DF[-2], 
                 FUN = function(X) paste(unique(X), collapse=", "))   
colnames(DF2) <- c("protein.acs", "Go.Mf")
compl <- merge(pm.5767, DF2, by="protein.acs")
grep("GPI", compl$Go.Mf)
write.csv2(compl, file="~/MP_RStudio_May2014/Aminooxybiotin/P5767_P5768/P5767/PM_P5767_Mf.csv")

# First tryout of Silica beads =====================================================================
# P5572
pm.5572 <- read.csv2("~/MP_RStudio_May2014/SilicaBeads/PM_P5572.csv", stringsAsFactors=F, header=T)
pmonly.5572 <- read.csv2("~/MP_RStudio_May2014/SilicaBeads/PMOnly_P5572.csv", stringsAsFactors=F, header=T)

length(grep("integral", pmonly.5572$Go.Cc))
to.assign <- setdiff(pm.5572$Protein.Acs, pmonly.5572$Protein.Acs)
new.df <- pm.5572[grep(paste(to.assign, collapse="|"), pm.5572$Protein.Acs), ]
length(grep("extracellular", new.df$Go.Cc))
length(grep("cell surface", new.df$Go.Cc))
length(grep("glycoprotein", new.df$Description))

go.terms=getBM(attributes=c("accession", "db2go_f__dm_description"),filters="accession",
               values=pm.5572$Protein.Acs,mart=unimart)

# Organize the output as dataframe
DF <- as.data.frame(go.terms)

# Improve the representation to show one protein per row
DF2 <- aggregate(DF[2], DF[-2], 
                 FUN = function(X) paste(unique(X), collapse=", "))   
colnames(DF2) <- c("Protein.Acs", "Go.Mf")
compl <- merge(pm.5572, DF2, by="Protein.Acs")
grep("GPI", compl$Go.Mf)
write.csv2(compl, file="~/MP_RStudio_May2014/SilicaBeads/PM_P5572_Mf.csv")

# Second tryout Silica Beads ============================================================================
# P5745
pm.5745 <- read.csv2("~/MP_RStudio_May2014/SilicaBeads/PM_P5745.csv", stringsAsFactors=F, header=T)
pmonly.5745 <- read.csv2("~/MP_RStudio_May2014/SilicaBeads/PMOnly_P5745.csv", stringsAsFactors=F, header=T)

length(grep("integral", pmonly.5745$Go.Cc))
to.assign <- setdiff(pm.5745$protein.acs, pmonly.5745$protein.acs)
new.df <- pm.5745[grep(paste(to.assign, collapse="|"), pm.5745$protein.acs), ]
length(grep("extracellular", new.df$Go.Cc))
length(grep("cell surface", new.df$Go.Cc))
length(grep("glycoprotein", new.df$Description))

go.terms=getBM(attributes=c("accession", "db2go_f__dm_description"),filters="accession",
               values=pm.5745$protein.acs,mart=unimart)

# Organize the output as dataframe
DF <- as.data.frame(go.terms)

# Improve the representation to show one protein per row
DF2 <- aggregate(DF[2], DF[-2], 
                 FUN = function(X) paste(unique(X), collapse=", "))   
colnames(DF2) <- c("protein.acs", "Go.Mf")
compl <- merge(pm.5745, DF2, by="protein.acs")
grep("GPI", compl$Go.Mf)
write.csv2(compl, file="~/MP_RStudio_May2014/SilicaBeads/PM_P5745_Mf.csv")

# Biotin elution =======================================================================================
# P6157
pm.6157 <- read.csv2("~/MP_RStudio_May2014/CellSurfaceExtractionKit/P6157/PM_6157.csv", stringsAsFactors=F, header=T)
pmonly.6157 <- read.csv2("~/MP_RStudio_May2014/CellSurfaceExtractionKit/P6157/PMOnly_P6157.csv", stringsAsFactors=F, header=T)

length(grep("integral", pmonly.6157$Go.Cc))
to.assign <- setdiff(pm.6157$protein.acs, pmonly.6157$protein.acs)
new.df <- pm.6157[grep(paste(to.assign, collapse="|"), pm.6157$protein.acs), ]
length(grep("extracellular", new.df$Go.Cc))
length(grep("cell surface", new.df$Go.Cc))
length(grep("glycoprotein", new.df$Description))

go.terms=getBM(attributes=c("accession", "db2go_f__dm_description"),filters="accession",
               values=pm.6157$protein.acs,mart=unimart)

# Organize the output as dataframe
DF <- as.data.frame(go.terms)

# Improve the representation to show one protein per row
DF2 <- aggregate(DF[2], DF[-2], 
                 FUN = function(X) paste(unique(X), collapse=", "))   
colnames(DF2) <- c("protein.acs", "Go.Mf")
compl <- merge(pm.6157, DF2, by="protein.acs")
grep("GPI", compl$Go.Mf)
write.csv2(compl, file="~/MP_RStudio_May2014/CellSurfaceExtractionKit/P6157/PM_6157_Mf.csv")

# P5832 ================================================================================================
pm.5832 <- read.csv2("~/MP_RStudio_May2014/CellSurfaceExtractionKit/P5832/PM_5832.csv", stringsAsFactors=F, header=T)
pmonly.5832 <- read.csv2("~/MP_RStudio_May2014/CellSurfaceExtractionKit/P5832/PMOnly_P5832.csv", stringsAsFactors=F, header=T)

length(grep("integral", pmonly.5832$Go.Cc))
to.assign <- setdiff(pm.5832$protein.acs, pmonly.5832$protein.acs)
new.df <- pm.5832[grep(paste(to.assign, collapse="|"), pm.5832$protein.acs), ]
length(grep("extracellular", new.df$Go.Cc))
length(grep("cell surface", new.df$Go.Cc))
length(grep("glycoprotein", new.df$Description))

go.terms=getBM(attributes=c("accession", "db2go_f__dm_description"),filters="accession",
               values=pm.5832$protein.acs,mart=unimart)

# Organize the output as dataframe
DF <- as.data.frame(go.terms)

# Improve the representation to show one protein per row
DF2 <- aggregate(DF[2], DF[-2], 
                 FUN = function(X) paste(unique(X), collapse=", "))   
colnames(DF2) <- c("protein.acs", "Go.Mf")
compl <- merge(pm.5832, DF2, by="protein.acs")
grep("GPI", compl$Go.Mf)
write.csv2(compl, file="~/MP_RStudio_May2014/CellSurfaceExtractionKit/P6157/PM_6157_Mf.csv")

# P6287 ===============================================================================================
pm.6287 <- read.csv2("~/MP_RStudio_May2014/CellSurfaceExtractionKit/P6287/PM_P6287.csv", stringsAsFactors=F, header=T)
pmonly.6287 <- read.csv2("~/MP_RStudio_May2014/CellSurfaceExtractionKit/P6287/PMOnly_P6287.csv", stringsAsFactors=F, header=T)

length(grep("integral", pmonly.6287$Go.Cc))
to.assign <- setdiff(pm.6287$protein.acs, pmonly.6287$protein.acs)
new.df <- pm.6287[grep(paste(to.assign, collapse="|"), pm.6287$protein.acs), ]
length(grep("extracellular", new.df$Go.Cc))
length(grep("cell surface", new.df$Go.Cc))
length(grep("glycoprotein", new.df$Description))

go.terms=getBM(attributes=c("accession", "db2go_f__dm_description"),filters="accession",
               values=pm.6287$protein.acs,mart=unimart)

# Organize the output as dataframe
DF <- as.data.frame(go.terms)

# Improve the representation to show one protein per row
DF2 <- aggregate(DF[2], DF[-2], 
                 FUN = function(X) paste(unique(X), collapse=", "))   
colnames(DF2) <- c("protein.acs", "Go.Mf")
compl <- merge(pm.6287, DF2, by="protein.acs")
grep("GPI", compl$Go.Mf)
write.csv2(compl, file="~/MP_RStudio_May2014/CellSurfaceExtractionKit/P6287/PM_P6287_Mf.csv")

=========================================================================================================
# Make a stacked barplot for that
totals <- c(881,477,440,394,394,263,274,258,129)
PMOnly <- c(471,222,213,190,196,178,90,91,48)
int <- c(264,113,115,103,100,133,23,17,14)
ext <- c(410,255,227,204,198,85,184,167,81)
other.pm <- PMOnly-int

vec <- rbind(int, other.pm, ext)
colnames(vec) <- c("Biotin-SDS 1", "Biotin-SDS 2", "Biotin 1", "Biotin 2", "Biotin 3", "Aminooxy 1", "Aminooxy 2", "Silica 1", "Silica 2")
mat <- as.matrix(vec)

pdf("~/MP_RStudio_May2014/Triton_X-114/StackedBarplot_PhaseDistrib.pdf")
barplot(mat, col=cm.colors(2), cex.names = 1, las=1, xlim=c(0, 10.5), ylab="# of proteins")
legend("bottomleft", legend = c("Membrane", "Other"), fill = cm.colors(2))
dev.off();

require(reshape2)
melt.mat <- melt(mat)

require(ggplot2)
c <- ggplot(melt.mat, aes(x=Var2, y=value, fill=Var1))
c+geom_bar(stat="identity")+ylab("# of PM proteins")+xlab("")+scale_fill_brewer(palette="Paired")
