library("NetRep")

#import data
setwd("/media/jack/SP PHD U3/Network analysis/Network Analysis PD (newest)/1 Sorting data to be used")
load("HCdatExpr.Rdata")
load("PDdatExpr.Rdata")



####correlation in test and discovery datasets
discovery_data <- HCdatExpr
test_data <- PDdatExpr
data_list <- list(cohort1=discovery_data, cohort2=test_data)

discovery_correlation <- cor(data_list$cohort1)
test_correlation <- cor(data_list$cohort2)
correlation_list <- list(cohort1=discovery_correlation, cohort2=test_correlation)



####import adjacencies
setwd("/media/jack/SP PHD U3/Network analysis/Network Analysis PD (newest)/2 WGCNA/HC")
load("adjacency.Rdata")
discovery_network <- adjacency #(from HC)
rm(adjacency)

setwd("/media/jack/SP PHD U3/Network analysis/Network Analysis PD (newest)/2 WGCNA/PD")
load("adjacency.Rdata")
test_network <- adjacency #(from PD)
network_list <- list(cohort1=discovery_network, cohort2=test_network)
rm(adjacency)



#import the modules
setwd("/media/jack/SP PHD U3/Network analysis/Network Analysis PD (newest)/3 K-means correction of data/HC")
load("adjusted_module_info.Rdata")
module_labels <- moduleColors 



#calculate module preservation
preservation <- modulePreservation( network=network_list, data=data_list, correlation=correlation_list, moduleAssignments=module_labels, discovery="cohort1", test="cohort2", nPerm=10000, nThreads=5, alternative = "less" )



setwd("/media/jack/SP PHD U3/Network analysis/Network Analysis PD (newest)/4 Identify changed modules/Corrected/Using Lesser/HC not preserved in PD")
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
rownames(ModulePreservation) <- rownames(preservation$p.value)

ModulePreservation[,1] <- apply(preservation$p.value, 1, min)
ModulePreservation[,2] <- apply(FDR.p.values, 1, min)
ModulePreservation[,3] <- apply(Bonf.p.values, 1, min)

ModulePreservation <- ModulePreservation[order(ModulePreservation[,3], ModulePreservation[,2], ModulePreservation[,1],  decreasing = FALSE),]

save(ModulePreservation, file = "ModulePreservation.Rdata")
