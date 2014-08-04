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
	for( i in 1:( ncol(data)-1))
	{
		Tuple <- aData$getMaximumSimilarity()
		#append(tree$genericPairs, aData$findGenericPairAtLevel(Tuple))
		#tree$appendHierarchyLevel(Tuple)
		aData$joinClasses(Tuple)
	}

	if(report){
		report = list()
		report$initialSimilarityMatrix<-aData$getSimilarityMatrixAtLevel(0)
		report$cor <-cor(data)
		report$basicStats <- cbind(apply(data,2,sum), apply(data, 2, mean), apply(data,2, sd))
		colnames(report$basicStats)<-c("freq","mean","sd")

	}
	#tree$significative <- ASImodel$getSignificativeNodes()
	#return(tree)
}
