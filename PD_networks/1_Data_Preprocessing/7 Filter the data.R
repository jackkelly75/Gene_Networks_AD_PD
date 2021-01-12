########to sort
library(affy)
library(hgu133plus2.db)


########################################################
#filter out any probes with duplicate genes or no gene #
########################################################

#all data
probes=row.names(data.rma.ComBat)     #get a list of probes
Symbols = unlist(mget(probes, hgu133plus2SYMBOL, ifnotfound=NA))   #get a list of symbols that match the probes of raw.data
Entrez_IDs = unlist(mget(probes, hgu133plus2ENTREZID, ifnotfound=NA))   #get a list of entrez IDs that match the probes of raw.data
data.rma.ComBat.annotate =cbind(probes,Symbols,Entrez_IDs,data.rma.ComBat) #binds probes, entrez Ids, symbols and expression set together

remove = function(x) {
  comma = grep(",", rownames(x))    #find any row with comma in them
        if(length(comma)>0){
          x = x[-comma,]
        }
  space = grep(" ", rownames(x))    #find any row with space in them
        if(length(space)>0){
          x = x[-space,]
        }
  semicolon = grep(";", rownames(x))    #find any row with semicolon in them
        if(length(semicolon)>0){
          x = x[-semicolon,]
  }
  colon = grep(":", rownames(x))    #find any row with colon in them
        if(length(colon)>0){
          x = x[-colon,]
        }
  i <- is.na(x[,3])
  l <- x[!i,]      # remove any probes that don't have a mapped gene
  id=grep("AFFX",rownames(l))
  m = l[-id,]      #remove any probes that are controlr AFFX probes
}

data.rma.ComBat.sorted = remove(data.rma.ComBat.annotate)
data.rma.ComBat.sorted = data.rma.ComBat[rownames(data.rma.ComBat.sorted),]
nrow(data.rma.ComBat.sorted) # 41925



################################################################
#If multiple probes map to one gene, only keep highest MAD one #
################################################################
#All data has 434 samples and 41925 probes

#All data
MAD <- vector(mode="numeric", length=0)
for( i in 1:41925){                
  MAD[i] <- mad(data.rma.ComBat.sorted[i,1:434])
}
ExprsMAD <- cbind(MAD, data.rma.ComBat.sorted)


#convert to data frame to allow different data types in one table
ExprsMAD <- as.data.frame(ExprsMAD)

##Annotating
#All data
probes=row.names(ExprsMAD)
Entrez_IDs = unlist(mget(probes, hgu133plus2ENTREZID, ifnotfound=NA)) 
AllData = cbind(Entrez_IDs, ExprsMAD)

##can then use code like this to select the genes with highest mad if are entrez same
uniqGeneId <- function(x) {
   x = x[order(x[,1], abs(x[,2]), decreasing = TRUE), ]
   entrezID = unique(x[,1])
   id = match(entrezID, x[,1])
   x = x[id[!is.na(id)], ]
   x
}

dim(AllData[duplicated(AllData[,1]),])[1]
All_data <- uniqGeneId(AllData)   # 20186 probes
All_data <- All_data[,3:ncol(All_data)]


#########################################################
#filtering by mean - with the bottom 5% being removed   #
#########################################################
#for all the data
myVar <- rowMeans(All_data)
myVar1 <- myVar[myVar > quantile(myVar, 0.05) ]  
filtered.data = All_data[names(myVar1),]

nrow(filtered.data) #19176   #all this length

######################################
#annotate the data                   #
######################################
#Set row names as gene Symbols
#All data
probes=row.names(filtered.data)
rownames(filtered.data) = unlist(mget(probes, hgu133plus2SYMBOL, ifnotfound=NA)) 

##########################################
#transpose the data to be used in WGCNA  #
##########################################
datExpr <- t(filtered.data)
