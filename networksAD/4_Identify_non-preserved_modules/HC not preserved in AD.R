library(NetRep)


#import data
setwd("D:/Network analysis/Network Analysis AD/GSE63060 + GSE63061/1 Sorting data to be used")
load("CTLdatExpr.Rdata")



####correlation in test and discovery datasets
discovery_data <- CTLdatExpr
load("ADdatExpr.Rdata")
test_data <- ADdatExpr
data_list <- list(cohort1=discovery_data, cohort2=test_data)

discovery_correlation <- cor(data_list$cohort1)
test_correlation <- cor(data_list$cohort2)
correlation_list <- list(cohort1=discovery_correlation, cohort2=test_correlation)



####import adjacencies
setwd("D:/Network analysis/Network Analysis AD/GSE63060 + GSE63061/2 WGCNA/HC")
load("adjacency.Rdata")
discovery_network <- adjacency #(from CTL)
rm(adjacency)
setwd("D:/Network analysis/Network Analysis AD/GSE63060 + GSE63061/2 WGCNA/AD")
load("adjacency.Rdata")
test_network <- adjacency #(from AD)
rm(adjacency)
network_list <- list(cohort1=discovery_network, cohort2=test_network)
rm(discovery_network)
rm(test_network)



#import the modules
setwd("D:/Network analysis/Network Analysis AD/GSE63060 + GSE63061/3 K-means correction of data/HC")
load("adjusted_module_info.Rdata")
module_labels <- moduleColors ##from the AD merged



#calculate module preservation
preservation <- modulePreservation( network=network_list, data=data_list, correlation=correlation_list, moduleAssignments=module_labels, discovery="cohort1", test="cohort2", nPerm=10000, nThreads=5, alternative = "less")

setwd("D:/Network analysis/Network Analysis AD/GSE63060 + GSE63061/4 Identify changed modules/Corrected/Lesser")
save(preservation, file = "preservation.Rdata")


######################################
#organise into corrected and table   #
######################################

FDR.p.values <- preservation$p.values
Bonf.p.values <- preservation$p.values

for (i in 1:7){
	FDR.p.values[,i] <- p.adjust(FDR.p.values[,i], method = "fdr")
	Bonf.p.values[,i] <- p.adjust(FDR.p.values[,i], method = "bonferroni")
}


ModulePreservation <- matrix(, nrow = length(unique(module_labels)), ncol = 3)

colnames(ModulePreservation) <- c("Highest pvalue", "FDR corrected", "Bonferroni corrected")

df <- data.frame(unique(moduleLabels),unique(moduleColors))
df[order(df[,1]),]
df <- df[order(df[,1]),]
rownames(preservation$p.value) <- df[,2]
rownames(ModulePreservation) <- df[,2]

ModulePreservation[,1] <- apply(preservation$p.value, 1, min)
ModulePreservation[,2] <- apply(FDR.p.values, 1, min)
ModulePreservation[,3] <- apply(Bonf.p.values, 1, min)

ModulePreservation <- ModulePreservation[order(ModulePreservation[,3], ModulePreservation[,2], ModulePreservation[,1],  decreasing = FALSE),]

save(ModulePreservation, file = "ModulePreservation.Rdata")
