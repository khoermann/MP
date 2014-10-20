# Making a stacked bar plot

# For the Biotinylation "efficiency", thus how many proteins carried a detected biotinylation modific. 
totals <- c(966,821,729,671,646)
biotinylated <- c(162,81,152,85,38)
non.biotinylated <- totals-biotinylated
vec <- rbind(biotinylated, non.biotinylated)
colnames(vec) <- c("classic", "1st tryout", "2nd tryout", "3rd tryout", "Urea wash")
mat <- as.matrix(vec)

library(RColorBrewer)
sequential <- brewer.pal(2, "YlOrRd")
barplot(mat, col=sequential, cex.names = 1, las=1, xlim=c(0, 7.5), ylab="# of proteins")
legend("bottomright", legend = c("YES", "NO"), fill = sequential[1:2], title='Biotinylation status')


# Now only look at all PM proteins and see what fraction thereof is biotinylated
PM <- c(477,440,394,394,380)
PM.biotinylated <- c(54,38,62,41,17)
non.biotinylated <- PM-PM.biotinylated
vec <- rbind(PM.biotinylated, non.biotinylated)
colnames(vec) <- c("classic", "1st tryout", "2nd tryout", "3rd tryout", "Urea wash")
mat <- as.matrix(vec)

library(RColorBrewer)
sequential <- brewer.pal(2, "YlOrRd")
barplot(mat, col=sequential, cex.names = 1, las=1, xlim=c(0, 7.5), ylab="# of PM proteins")
legend("bottomright", legend = c("YES","NO"), fill = sequential[1:2], title='Biotinylation status')



