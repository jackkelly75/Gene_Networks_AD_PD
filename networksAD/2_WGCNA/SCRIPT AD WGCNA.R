library(WGCNA)
enableWGCNAThreads()
library(flashClust)


##########
# Identify the best soft thresholding power
##########
# Create table of soft threshold powers against scale free network topology (signed R2)
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(ADdatExpr, powerVector = powers, verbose = 5, networkType = "signed")
#indicates a softpower of 4 as the best

# Plot the soft threshold powers against scale free network topology (signed R2)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of 0.85
abline(h=0.85,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
##dev.off()


##############################################
#use the soft power of 4 to create network  #
##############################################

softPower = 4
##create adjacency matrix
adjacency = adjacency(ADdatExpr, power = softPower, type = "signed")

#######plot to see the scale free topology
k=softConnectivity(datE=ADdatExpr,power=softPower)
sizeGrWindow(10,5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k, main="Check scale free topology\n")

########create a topological overlap matrix
TOM = TOMsimilarity(adjacency, TOMType="signed")

##########turn into distance matrix
dissTOM = 1-TOM

geneTree = flashClust(as.dist(dissTOM), method="average") #flashClust is faster hclust

plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)

save(adjacency, file = "adjacency.Rdata")
save(TOM, dissTOM, geneTree, file = "TOMs and geneTree.Rdata")


###################################
#identify communities in network  #
###################################

minModuleSize = 10
#Adaptive Branch Pruning Of Hierarchical Clustering Dendrograms
dynamicMods = cutreeDynamic(dendro= geneTree, distM= dissTOM, deepSplit=3, pamRespectsDendro= FALSE, minClusterSize = minModuleSize)
# convert the labels to colours
dynamicColors= labels2colors(dynamicMods)

# plot the dendrograms
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")

# Calculates module eigengenes of modules
MEList= moduleEigengenes(ADdatExpr, colors= dynamicColors,softPower = softPower)
MEs= MEList$eigengenes
MEDiss= 1-cor(MEs)
METree= flashClust(as.dist(MEDiss), method= "average")

# plot the clustering of module eigengenes
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
MEDissThres = 0.2
abline(h=MEDissThres, col = "blue")
MEDissThres = 0.05
abline(h=MEDissThres, col = "red")

save(dynamicMods, dynamicColors, MEList, MEs, file = "InitialCommunities.Rdata")



########################################
#merge any close modules               #
########################################
MEDissThres = 0.05
merge = mergeCloseModules(ADdatExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
#drops 31 to 29 modules

# merged module colors
moduleColors = merge$colors
names(moduleColors) <- colnames(ADdatExpr)
#eigengenes of the new merged modules:
mergedMEs = merge$newMEs


# create list of colours starting with grey and then assign colours
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
names(moduleLabels) <- colnames(ADdatExpr)
MEs=mergedMEs
rownames(MEs)<-rownames(ADdatExpr)

save(MEs, moduleLabels, moduleColors , file = "AfterMergedMEs.Rdata")

#now we want to calculate the intramodular connectivity, how well specific genes are connected to other genes within same module
Alldegrees1=intramodularConnectivity(adjacency, moduleColors)
save(Alldegrees1, file = "Alldegrees1.Rdata")


### plot dendrograms
plotDendroAndColors(geneTree, cbind(dynamicColors, moduleColors), c("Dynamic Tree Cut", "Merged dynamic"),dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)



################################
#print the modules out to csv  #
################################
colors <- unique(moduleColors)
for (i in 1:length(colors)){
	MyData <- names(which(moduleColors == colors[i]))
	write.csv(MyData, file = paste(colors[i],".csv", sep = ""))
}

Gene_numbers <- vector(mode="numeric", length=0)
for (i in 1:length(colors)){
	Gene_numbers[i] <- length(which(moduleColors == colors[i]))
}
names(Gene_numbers) <- colors
save(Gene_numbers, file = "Gene_numbers.Rdata")
