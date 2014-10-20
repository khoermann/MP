## Deeper analysis of most successful tryout of Silica beads protocol (30 mio. cells)

  getwd(); setwd('/slow/home/khoermann/MP/SilicaBeads/')
  library(xlsx)

P5571.df <- read.xlsx2('20131203.SilicaBeads.P5571.results.xlsx', sheetIndex=1)


# Search results with a vector of GOCC annotations of interest
toMatch <- c('plasma membrane', 'cell surface', 'extracellular')
matches.pm <- P5571.df[grep(paste(toMatch, collapse='|'), P5571.df$Go.Cc), ]
  
# Identify non-PM proteins  
matches.not.pm <- P5571.df[grep(paste(toMatch, collapse='|'), P5571.df$Go.Cc, invert=TRUE), ]

# Extract Global Spectra Counts for the two groups 
  # Turn Global Spectra Counts into a numeric vector
  gsc.pm <- as.numeric(as.character(matches.pm$Global.Spectra.Count))
  gsc.not.pm <- as.numeric(as.character(matches.not.pm$Global.Spectra.Count))
  
  # Plot the data
  par(mfrow=c(1,1))
  plot(gsc.not.pm, cex=0.75, ylab='Global Spectra Count')
  points(gsc.pm, col='red', cex=0.75)
  legend(locator(1), c('Non-PM', "PM"), fill=c('black', 'red'), cex=rep(0.75,2))
  
# Display data distribution graphically
  par(mfrow=c(1,1))
  boxplot(gsc.pm, gsc.not.pm, ylab='Global Spectra Count', names=c('PM', 'Non-PM'))

par(mfrow=c(1,2))  
hist(gsc.pm,
     ylim=c(0,100),
     xlab='Global Spectra Count',
     main='Histogram of PM protein distribution')
hist(gsc.not.pm,
     xlab='Global Spectra Count',
     main='Histogram of non-PM protein distribution')
   
# Get parameters for the two groups
summary(gsc.pm)
summary(gsc.not.pm)
