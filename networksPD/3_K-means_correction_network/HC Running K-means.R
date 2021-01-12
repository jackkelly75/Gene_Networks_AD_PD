# additional k-means step based on following paper:
#https://bmcsystbiol.biomedcentral.com/articles/10.1186/s12918-017-0420-6
library(km2gcn)
library(WGCNA)

# import data
setwd("D:/Network analysis/Network Analysis PD (newest)/2 WGCNA/HC")
load("AfterMergedMEs.Rdata")
load("TOMs and geneTree.Rdata")  #to get TOM
setwd("D:/Network analysis/Network Analysis PD (newest)/1 Sorting data to be used")
load("HCdatExpr.Rdata")



#Step 1. Let D be the expression data in which dij in D represents the expression value for sample i and gene j, being s samples and g genes in total.
 
net.file = list(moduleColors, MEs)
names(net.file) <- c("moduleColors", "MEs")

expr.data = HCdatExpr
beta=10
n.iterations=400
meg = 0
tom.matrix=TOM
plot.evolution=TRUE
plot.go=FALSE
debug=NULL
net.type="signed"
min.genes.for.grey=20

cat("Working with",nrow(expr.data),"samples and",ncol(expr.data),"genes\n")
##Working with 230 samples and 19176 genes




#Step 2. Construct the partition by the WGCNA process, let P_D={m_1, m_2, ..., m_n} be that partition where m_k is the k-th module.

net = net.file
cat("The network includes",length(net$moduleColors),"genes and ", length(unique(net$moduleColors)),"modules\n")
#The network includes 19176 genes and  25 modules

partition.in.colors <- net$moduleColors



#Step 3. Get the eigengenes for each module within the partition, E={e_1, e_2, ..., e_n}

if(sum(partition.in.colors == "grey") < min.genes.for.grey){
	eigengenes = moduleEigengenes(expr.data,partition.in.colors, excludeGrey=TRUE)
}else{
	eigengenes = moduleEigengenes(expr.data,partition.in.colors, excludeGrey=F)
}

cat("We got",length(eigengenes$eigengenes)," eigengene vectors\n")
#We got 25  eigengene vectors

centroid.labels <- substring(names(eigengenes$eigengenes),3)
print("Module colors are")
print(centroid.labels)





#Step 4. Set up the k-means clustering
#Step 4.1. Set k to n
#Step 4.2. Set the centroids C to the eigengenes E, thus C to E

k <- length(eigengenes$eigengenes)

createCentroidMatrix= function(centroids, features)
{
  centroidMatrix = matrix(data=0,nrow = features, ncol=length(centroids))
  for(i in 1:k)
  {
    centroidMatrix[,i]= dataMatrix[, centroids[i]]
  }
  return(centroidMatrix)
}

centroids <- matrix(unlist(eigengenes$eigengenes), ncol = ncol(eigengenes$eigengenes))
colnames(centroids) <- colnames(eigengenes$eigengenes)




#Step 5. Run the algorithm and monitor its evolution
#Step 5.1 Set iterations to 0
#Step 5.2 Create a new partition P', given C with n modules such that, for each gene, 1 <= j <= g, g_j belongs to the module c_t in C such that a distance meassure d(g_j,c_t) is minimum.
#Step 5.3 Calculate eigengenes of P', giving a new E'
#Step 5.4 Evaluate the progress. If progress done, set iterations to iterations + 1 and C to E' and go to step 5.2
#Step 5.5 Finish

partitions <- list()
new.partition <- match(partition.in.colors, centroid.labels)
names(new.partition) <- centroid.labels[new.partition]
partitions[[1]] <- new.partition


#Launch the iterations
exchanged.genes = meg + 1
iteration = 1


getBestModuleCor <- function(gene,centroids,signed=TRUE){
  return(which.max(0.5 * (1 + cor(centroids, gene, use = "all.obs", method = "pearson"))))
}

getExchangedGenes <- function(old.partition,new.partition){   		stopifnot(length(old.partition) == length(new.partition))
	return(old.partition[old.partition != new.partition])
}

getNewCentroids <- function(expr.data,partition.in.colors,centroid.labels,mgg){
  if(sum(partition.in.colors == "grey") < mgg)
    eg.vectors = moduleEigengenes(expr.data,partition.in.colors, excludeGrey=TRUE)$eigengenes
  else
    eg.vectors = moduleEigengenes(expr.data,partition.in.colors, excludeGrey=F)$eigengenes

  names(eg.vectors) <- substring(names(eg.vectors),3)
  eg.vectors <- eg.vectors[,centroid.labels]
  return(eg.vectors)
}


while(exchanged.genes > meg & iteration <= n.iterations){
    print(paste0("Starting partition ",iteration))
    print(paste0("Number of centroids before getting new partition ",ncol(centroids)))
    new.partition <- apply(expr.data,MARGIN=2,getBestModuleCor,centroids=centroids,signed=(net.type == "signed"))
    partitions[[iteration + 1]] <- new.partition
    exchanged.gene.count <- length(getExchangedGenes(partitions[[iteration]], partitions[[iteration + 1]]))
    cat("A total of ", exchanged.gene.count, " genes moved to another partition\n")
    new.partition.in.colors <- centroid.labels[unlist(new.partition)]
    centroids <- getNewCentroids(expr.data,new.partition.in.colors,centroid.labels,min.genes.for.grey)
    exchanged.genes = exchanged.gene.count
    iteration = iteration + 1
}

cat("We finish with",iteration,"iterations\n")
#71
cat("Last number of gene changes where",exchanged.genes,"\n")
save(partitions, file = "partitions.Rdata")


moduleLabels <- partitions[[71]]
names(new.partition.in.colors) <- names(moduleLabels)
moduleColors <- new.partition.in.colors
MEs <- centroids

save(moduleLabels, moduleColors, MEs, file = "adjusted_module_info.Rdata")



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



#############################################
#plot graph of genes moved each iteration   #
#############################################
net.label = "HC Data genes Exhanged in k-means adjustment"
g.changes = NULL
n = length(partitions)
for(index in 1:(n-1)){
  g.changes = c(g.changes,sum(partitions[[index]] != partitions[[index + 1]]))
}
plot(1:(length(partitions)-1),g.changes, main=net.label ,xlab="Iteration number", ylab = "Number of genes moved")
