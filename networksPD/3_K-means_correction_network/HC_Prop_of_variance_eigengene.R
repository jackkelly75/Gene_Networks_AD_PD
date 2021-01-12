library(WGCNA)


expr <- HCdatExpr

colors <- moduleColors

temp <- moduleEigengenes(expr,colors,nPC = 10,excludeGrey = FALSE,softPower = 11)

HCvarExplain <- temp$varExplained
colnames(HCvarExplain) <- substr(colnames(temp$averageExpr), 3, 100)

save(HCvarExplain, file = "HCvarExplain.Rdata")
