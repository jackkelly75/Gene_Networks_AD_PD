#########################################################
#set up the phenodata files after removing the outliers #
#########################################################
#to remove
#GSM2630908 - PINK1 mutation
#GSM2630938
#GSM2631142
#GSM2630912

remove = c("GSM2630908", "GSM2630938", "GSM2631142", "GSM2630912")
pData <- pData[!rownames(pData) %in% remove, ]

HCpData <- pData[pData$`disease label:ch1`  == "CONTROL",]
#nrow(HCpData)  #230
PDpData <- pData[pData$`disease label:ch1`  == "IPD",]
#nrow(PDpData)  #204
pData <- rbind(HCpData, PDpData)


################################
#read in the CEL data files    #
################################
cels = pData$fileName
unfiltered.raw.data = ReadAffy(filenames = cels, phenoData = pData)
