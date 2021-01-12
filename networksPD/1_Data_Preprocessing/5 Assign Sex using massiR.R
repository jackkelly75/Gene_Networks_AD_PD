#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4080740/

#############################
#Get missing sex data    #
#############################
library(massiR)
hgu133Plus2.probes <- data.frame(y.probes["affy_hg_u133_plus_2"])

# use the massi.y function to calculate probe variation
massi_y_out <- massi_y(expression_data= data.rma.norm, y_probes=hgu133Plus2.probes)

# plot probe variation to aid in deciding on the most informative subset of y chromosome probes
massi_y_plot(massi_y_out)

# Extract the informative probes for clustering
massi_select_out <- massi_select(data.rma.norm, hgu133Plus2.probes, threshold=4)

# cluster samples to predict the sex for each sample
massi_cluster_out <- massi_cluster(massi_select_out)
# get the predicted sex for each sample
gender_prediction <- data.frame(massi_cluster_out[[2]])
rownames(gender_prediction) <- gender_prediction[,1]
gender_prediction <- gender_prediction[,2:ncol(gender_prediction)]

##########################
#testing the results     #
##########################
trial <-merge(pData, gender_prediction, by = "row.names")

sextest <- trial[,c("sex", "Sex:ch1")]

sextest$sex <- gsub("female", "F", sextest$sex)
sextest$sex <- gsub("male", "M", sextest$sex)
sextest$`Sex:ch1` <- gsub("Male", "M", sextest$`Sex:ch1`)
sextest$`Sex:ch1` <- gsub("Female", "F", sextest$`Sex:ch1`)

sextest <- na.omit(sextest)
sextest$hey <- NA
for (i in 1:nrow(sextest)){
	sextest$hey[i] <- paste(sextest$`Sex:ch1`[i], sextest$sex[i], sep = "")
}
sextest <- sextest[sextest$hey != "MM",]
sextest <- sextest[sextest$hey != "FF",]
#sextest now contains the number of mismatching genders

#14 different between them, an error of 3.5%  ((14/399)*100)


##################################
#re creating pData with sex   #
##################################
pData <-merge(pData, gender_prediction, by = "row.names")
rownames(pData) <- pData[,1]
pData <- pData[,2:ncol(pData)]


pData$sex <- gsub("male", "Male", pData$sex)
pData$sex <- gsub("feMale", "Female", pData$sex)


x <- pData$`Sex:ch1`

y <- pData$sex

for (i in 1:length(y)){
	if(is.na(x[i])){
		x[i] <- y[i]
	}
}
pData$final_gender <- x
