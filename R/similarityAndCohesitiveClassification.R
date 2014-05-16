performASISimilarityClassification <-function(data, model = 'pois'){
	tree <- similarityTree$new()
	return(callASIalgorithm(data, tree, model))
}

# performASICohesitiveClassification <-function(data, model='pois'){
# 	tree <- cohesionTree$new()
# 	return(callASIAlgorithm( data, tree, model))
# }


callASIAlgorithm <-function( data, tree, model){
	ASIdata <- ASIdata$new( data, model)
	
	for( i in 1:( ncol(data)-1))
	{
		Tuple <- ASIdata.getMaximumSimilarity()
		append(tree$genericPairs, ASIdata$findGenericPairAtLevel(Tuple))
		tree$appendHierarchyLevel(Tuple)
		ASIdata$joinClasses(Tuple)
		ASIdata$computeSimilarityMatrix()
	}
	tree$significative <- ASIdata$getSignificativeNodes()
	return(tree)
}
