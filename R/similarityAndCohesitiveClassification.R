performASISimilarityClassification <-function(data, model = 'pois'){
	tree <- similarityTree$new()
	return(callASIalgorithm(data, tree, model))
}

# performASICohesitiveClassification <-function(data, model='pois'){
# 	tree <- cohesionTree$new()
# 	return(callASIAlgorithm( data, tree, model))
# }


callASIAlgorithm <-function( data, tree, model){
	aData <- ASImodel$new( data, model)

	for( i in 1:( ncol(data)-1))
	{
		Tuple <- ASImodel$getMaximumSimilarity()
		append(tree$genericPairs, aData$findGenericPairAtLevel(Tuple))
		tree$appendHierarchyLevel(Tuple)
		aData$joinClasses(Tuple)
		aData$computeSimilarityMatrix()
	}
	tree$significative <- ASImodel$getSignificativeNodes()
	return(tree)
}
