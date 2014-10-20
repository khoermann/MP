## Deeper analysis of SilicaBeads data

setwd('/slow/home/khoermann/MP/SilicaBeads/'); getwd(); 
  

silica.df <- read.csv2('20140325.SilicaBeads.P5744.P5745-34.GOterms.csv', header=TRUE, stringsAsFactors=F)
head(silica.df)
str(silica.df)
  
  
# Search results with a vector of GOCC annotations of interest
toMatch <- c('plasma membrane', 'cell surface', 'extracellular')
matches.pm <- silica.df[grep(paste(toMatch, collapse='|'), silica.df$go_name), ]
  length(unique(matches.pm[['accession']]))
  
write.csv2(matches.pm, file='PM_P5744.P5745-34.csv')  
  
double.check <- DF[grep(paste(toMatch, collapse='|'), DF$go_name), ]
  length(unique(double.check[['accession']]))

# Additional GO terms that are good to know
matches.itm <- silica.df[grep("integral to membrane", silica.df$go_name), ]
matches.mem <- silica.df[grep("membrane", silica.df$go_name), ]
# Intersect between membrane and PM
length(intersect(matches.mem$accession, matches.pm$accession))
  
# Identify non-PM proteins  
matches.not.pm <- silica.df[grep(paste(toMatch, collapse='|'), silica.df$go_name, invert=TRUE), ]

