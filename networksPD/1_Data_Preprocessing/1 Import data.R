library(hgu133plus2.db)
library(limma)
library(R.utils)
library(tidyr)
library(RColorBrewer)
library(calibrate)
library(GEOquery)
library(affy)

#############################
#create the phenodata file  #
#############################
gse99039 <- getGEO("GSE99039", GSEMatrix = TRUE)
pData <- pData(phenoData(gse99039[[1]]))[,c(1,34,48:54, 56, 57, 60:64)]
pData$fileName <- NA
for (i in 1:nrow(pData)){
	pData$fileName[i] <- paste(rownames(pData)[i], "_", pData$description[i], sep = "")
	pData$fileName[i] <- gsub("\\(|\\)", "", pData$fileName[i])
	pData$fileName[i] <- gsub("Plus_2.CEL", "Plus_2_.CEL", pData$fileName[i])
}

#############################
#set up the phenodata files #
#############################
HCpData <- pData[pData$`disease label:ch1`  == "CONTROL",]
PDpData <- pData[pData$`disease label:ch1`  == "IPD",]
pData <- rbind(HCpData, PDpData)


################################
#read in the CEL data files    #
################################
#untar("GSE99039_RAW.tar")
#cels = list.files(pattern = "CEL")
#sapply(paste(cels, sep= "/"), gunzip)
cels = pData$fileName
unfiltered.raw.data = ReadAffy(filenames = cels, phenoData = pData)
cels = PDpData$fileName
