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
		Tuple <- aData$getMaximumSimilarity()
		#append(tree$genericPairs, aData$findGenericPairAtLevel(Tuple))
		Tuple[Tuple>ncol(data)] = tree[Tuple[Tuple>ncol(data)]-ncol(data)]
		tree[i]=paste("(", Tuple[1],",",Tuple[2] ,")", sep="")
		aData$joinClasses(Tuple)
	}
	tree = paste(tree, ';', sep="")
	if(report){
		report = list()
		report$initialSimilarityMatrix<-aData$getSimilarityMatrixAtLevel(0)
		report$cor <-cor(data)
		report$basicStats <- cbind(apply(data,2,sum), apply(data, 2, mean), apply(data,2, sd))
		colnames(report$basicStats)<-c("freq","mean","sd")
		report$tree <- tree
	}
	#tree$significative <- ASImodel$getSignificativeNodes()
	#return(tree)
}
