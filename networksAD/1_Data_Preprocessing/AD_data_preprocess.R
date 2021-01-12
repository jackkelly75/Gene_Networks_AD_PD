###################
#load in sorted data
####################
##download the merged data
#avaible at https://figshare.com/s/78839db30d17d3f75aca

###read in the expression data table created
temp <- read.table("gse63060_61.merged.exp")
###convert to a matrix
info <- as.matrix(temp)

#set the colnames and then remove that row
colnames(info) <- info[1,]
info <- info[2:nrow(info),]

#set the rownames and then remove that coloumn
rownames(info) <- info[,1]
info <- info[,2:ncol(info)]
mode(info) <- "numeric" 

#rename and save
merged_data <- info
save(merged_data, file = "merged_data.Rdata")


############
#read in the info files
#############
gse63060 <- read.table("Samples_gse63060.info", sep="\t", header = T)
gse63061 <- read.table("Samples_gse63061.info", , quote = "", sep="\t", header = T, fill=TRUE)


##################
#sort them to include western european and other caucasian
##################
unique(gse63060$ethnicity)
##Western European, Other Caucasian,  Unknown 
gse63060 <- rbind(gse63060[gse63060$ethnicity == "Western European" ,], gse63060[gse63060$ethnicity == "Other Caucasian" ,])
rownames(gse63060) <- gse63060[,1]
gse63060 <- gse63060[,2:ncol(gse63060)]


unique(gse63061$ethnicity)
# Western European , Other Caucasian, Any_Other_White_Background, Caribbean, Irish , British  , Indian , British_English  , British English, British_Welsh, Any_Other_Asian_Background, White_And_Asian, British_Scottish, British_Other_Background, Any_Other_Ethnic_Background, Asian , Any_Other_Black_Background , unkown but she's white and speaks english with a slight south african accent
gse63061 <- rbind(gse63061[gse63061$ethnicity == "Western European" ,], gse63061[gse63061$ethnicity == "Other Caucasian" ,])


nrow(gse63060)   # 324
nrow(gse63060[gse63060$status == "AD",]) #143
nrow(gse63060[gse63060$status == "MCI",]) #77
nrow(gse63060[gse63060$status == "CTL",]) #104
#same as his paper

nrow(gse63061) #245
nrow(gse63061[gse63061$status == "AD",]) #102
nrow(gse63061[gse63061$status == "MCI",]) #65
nrow(gse63061[gse63061$status == "CTL",]) #78

#569 samples in total

merged.data <- merged_data[,c(rownames(gse63060), rownames(gse63061))]
pData <- rbind(gse63060[,1:5], gse63061[,1:5])
#status still has levels from the borderline data previously
pData$status <- droplevels(pData$status)



############################
#correcting all data       #
############################
library(sva)
pheno = pData
edata = merged.data
mod = model.matrix(~as.factor(status), data=pheno)

# reference-batch version, with covariates
batch = pheno$gender
data.rma.ComBat = ComBat(dat=edata, batch=batch,mod=mod,  par.prior=TRUE, prior.plots=FALSE)
batch = pheno$age
sorted_data = ComBat(dat=data.rma.ComBat, batch=batch, mod=mod, par.prior=TRUE, prior.plots=FALSE)


##############
#annotate the data and remove any NA
##############
library(illuminaHumanv3.db)
probes=row.names(sorted_data)
Symbols = unlist(mget(probes, illuminaHumanv3SYMBOL, ifnotfound=NA))
Entrez_IDs = unlist(mget(probes, illuminaHumanv3ENTREZID, ifnotfound=NA))
sorted_data.annotate =cbind(Entrez_IDs,Symbols, sorted_data)

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
	i <- is.na(x[,1])
	x <- x[!i,]      # remove any probes that don't have a mapped gene	  
}

sorted_data_annotate = remove(sorted_data.annotate)
#20019 genes for 14380 unique genes
##22756 was in his
sorted_annotate = sorted_data_annotate[,3:ncol(sorted_data_annotate)]
rownames(sorted_annotate) = sorted_data_annotate[,2]
mode(sorted_annotate) <- "numeric" 



################################################################
#If multiple probes map to one gene, only keep highest MAD one #
################################################################
#All data
MAD <- vector(mode="numeric", length=0)
for( i in 1:20019){                
        MAD[i] <- mad(sorted_annotate[i,1:569])
}


ExprsMAD <- cbind(MAD, sorted_annotate)
#convert to data frame to allow different data types in one table

uniqGeneId <- function(x) {
   x = x[order(rownames(x), abs(x[,2]), decreasing = TRUE), ]
   entrezID = unique(rownames(x))
   id = match(entrezID, rownames(x))
   x = x[id[!is.na(id)], ]
   x
}

ExprsMAD <- uniqGeneId(ExprsMAD)   # 14380 probes
ExprsMAD <- ExprsMAD[,2:ncol(ExprsMAD)]


#########################################################
#filtering by mean - with the bottom 5% being removed   #
#########################################################
#for all the data
myVar <- rowMeans(ExprsMAD)
myVar1 <- myVar[myVar > quantile(myVar, 0.05) ]  
datExpr = ExprsMAD[names(myVar1),]

nrow(datExpr) #13661


ADnames <- c(rownames(gse63060[gse63060$status == "AD",]), rownames(gse63061[gse63061$status == "AD",]))
CTLnames <- c(rownames(gse63060[gse63060$status == "CTL",]), rownames(gse63061[gse63061$status == "CTL",]))
MCInames <- c(rownames(gse63060[gse63060$status == "MCI",]), rownames(gse63061[gse63061$status == "MCI",]))

ADdatExpr <- t(datExpr[,ADnames])   #245 samples
CTLdatExpr <- t(datExpr[,CTLnames])   #182 samples
MCIdatExpr <- t(datExpr[,MCInames])   #142 samples
datExpr <- t(datExpr)
