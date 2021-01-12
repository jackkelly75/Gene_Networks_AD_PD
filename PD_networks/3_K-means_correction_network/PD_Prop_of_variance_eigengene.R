
expr <- PDdatExpr

colors <- moduleColors

temp <- moduleEigengenes(expr,colors,nPC = 10,excludeGrey = FALSE,softPower = 13)

PDvarExplain <- temp$varExplained
colnames(PDvarExplain) <- substr(colnames(temp$averageExpr), 3, 100)

save(PDvarExplain, file = "PDvarExplain.Rdata")
