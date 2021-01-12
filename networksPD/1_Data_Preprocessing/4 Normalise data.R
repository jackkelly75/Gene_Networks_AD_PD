library(affy)
library(RColorBrewer)


#############################
#normalise the data         #
#############################

#all data
data.rma.norm = rma(unfiltered.raw.data)
data.rma.expr = exprs(data.rma.norm)

########################
#plot normalised data  #
########################
brewer.cols <- brewer.pal(9, "Set1")

png("RMA normalised log scale probe intensities.png", width = 1200, height = 800)
par(mar=c(5.1,5.1,4.1,2.1))
boxplot(data.rma.expr, col = brewer.cols, ylab = "RMA normalised log (base 2) scale probe intensities", xlab = "Array names", main = "RMA normalised log (base 2) scale probe intensities of each array", xaxt='n', cex.main=2, cex.lab=1.5)
dev.off()
nrow(data.rma.expr) #54675
