#############################
#plot prenormalised plots   #
#############################
brewer.cols <- brewer.pal(9, "Set1")
png("Unprocessed log scale probe intensities.png", height = 700, width = 1200)
par(mar=c(5.1,5.1,4.1,2.1))
boxplot(unfiltered.raw.data, col = brewer.cols, ylab = "Unprocessed log (base 2) scale Probe Intensities", xlab = "Array Names", main = "Unprocessed log (base 2) scale Probe Intensities of each Array", xaxt='n', cex.main=2, cex.lab=1.5)
dev.off()
#Contruct density plots
png("Unprocessed density plots.png", height = 800, width = 1200)
par(mar=c(5.1,5.1,4.1,2.1))
hist(unfiltered.raw.data, main = "Density plot of log(2) probe intensities", col = brewer.cols, lty=1, xlab = "Log (base 2) Intensities", lwd =3, cex.main=2, cex.lab=1.5)
dev.off()


########################################
#getting the names of sample outliers  #
########################################
hist(PD.raw.data[,60:65], main = "PD Density plot of log(2) probe intensities", col = brewer.cols, lty=1, xlab = "Log (base 2) Intensities", lwd =3, cex.main=2, cex.lab=1.5)
samp.leg.names <- colnames(PD.raw.data[,60:65])
legend(12, 1.5, cex = 1.5, legend = samp.leg.names, lty = 1, col = brewer.cols, lwd = 3)
###GSM2630938


hist(HC.raw.data[,60:66], main = "HC Density plot of log(2) probe intensities", col = brewer.cols, lty=1, xlab = "Log (base 2) Intensities", lwd =3, cex.main=2, cex.lab=1.5)
samp.leg.names <- colnames(HC.raw.data[,60:66])
legend(12, 1.5, cex = 1.5, legend = samp.leg.names, lty = 1, col = brewer.cols, lwd = 3)
#GSM2630912


hist(HC.raw.data[,166:168], main = "HC Density plot of log(2) probe intensities", col = brewer.cols, lty=1, xlab = "Log (base 2) Intensities", lwd =3, cex.main=2, cex.lab=1.5)
samp.leg.names <- colnames(HC.raw.data[,166:168])
legend(12, 1.5, cex = 1.5, legend = samp.leg.names, lty = 1, col = brewer.cols, lwd = 3)
#GSM2631142
