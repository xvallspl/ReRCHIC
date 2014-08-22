callASIAlgorithm <-function( data, model, report=FALSE){
	aData <- ASImodel$new( data, model)
	tree = array('', ncol(data)-1)
	for( i in 1:( ncol(data)-1))
	{
		S <- aData$getMaximumSimilarity()
		#append(tree$genericPairs, aData$findGenericPairAtLevel(S$nodes))
		aData$joinClasses(S$nodes)
		S$nodes[S$nodes>ncol(data)] = tree[S$nodes[S$nodes>ncol(data)]-ncol(data)]
		tree[i]=paste("(", S$nodes[1],",",S$nodes[2] ,"):", S$value,sep="")
	}
	tree = paste(tree, ';', sep="")
	if(report){
		report = list()
		report$dims <- c(ncol(data), nrow(data))
		report$initialSimilarityMatrix<-aData$getSimilarityMatrixAtLevel(0)
		report$cor <-cor(data)
		report$freq <- NULL
		report$basicStats <- cbind(apply(data,2,sum), apply(data, 2, mean), apply(data,2, sd))
		colnames(report$basicStats) <- c("freq","mean","sd")
		report$tree <- tree
	}
}

displayReport <- function(report)
{
	print(paste("ncol:", report$dims[1], ", nrow:", report$dims[2]))
	print(report$basicStats)
	print("Bivariant frequency:")
	print("Correlation coefficient:")
	print("Similarity indices:")
	for(i in 1:ncol(length(report$tree)))
	{
		print(paste("Classification at level", i))
	}
}