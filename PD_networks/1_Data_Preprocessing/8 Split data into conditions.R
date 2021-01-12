PDpData <- pData[pData$`disease label:ch1` == "IPD",]
HCpData <- pData[pData$`disease label:ch1` == "CONTROL",]


PDdatExpr <- datExpr[rownames(PDpData),]
HCdatExpr <- datExpr[rownames(HCpData),]
