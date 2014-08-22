performASISimilarityClassification <-function(data, model = 'pois'){
	tree <- similarityTree$new()
	return(callASIalgorithm(data, tree, model))
}

# performASICohesitiveClassification <-function(data, model='pois'){
# 	tree <- cohesionTree$new()
# 	return(callASIAlgorithm( data, tree, model))
# }


callASIAlgorithm <-function( data, tree, model, report=FALSE){
	aData <- ASImodel$new( data, model)
	treeLevels = array('', ncol(data)-1)
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
		colnames(report$basicStats)<-c("freq","mean","sd")
		report$tree <- tree
	}
}

displayReport <- function(report)
{
	print(paste("ncol:", reports$dims[1], ", nrow:", reports$dims[2]))
	print(report$basicStats)
	print("Frecuencia de parejas de variables:")
	print("Coeficiente de correlación:")
	print("Indices de similaridad:")
	for(i in 1:ncol(length(report$tree)))
	{
		print(paste("Clasificación en el nivel", i))
	}
	
}