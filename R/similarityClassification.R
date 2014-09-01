#' Functions that performs the algorithm for the Statistical implicative analysis.
#'
#' @name callASIAlgorithm
#'
#' @param data 		External data.
#' @param model 	Probability model to use. Binom or pois.
#' @param report 	If true, displays a report
#' 
#' @return report 		Report of the data: algorithm results, hierarchical tree, basic statistic analysis. 
#' @return newickTree	Hierarchical tree in Newick's format.
#' 
#' @author Xavier Valls \email{xaviervallspla@@gmail.com}
#' @export

callASIAlgorithm <-function( data, model, report=FALSE){
	aData <- ASImodel$new(data)
	tree<-array("",ncol(data)-1)
	for( i in 1:( ncol(data)-1))
	{
		S <- aData$getMaximumSimilarity()
		aData$joinClasses(S$nodes)
		#genericPair<- aData$findGenericPairAtLevel(S$nodes)
		#aData$setGenericImplications(genericPair$pos, i)
		S$nodes[S$nodes> aData$nPrimitiveClasses] = tree[S$nodes[S$nodes>aData$nPrimitiveClasses]-aData$nPrimitiveClasses]
		tree[i]=paste("(", S$nodes[1],",",S$nodes[2] ,"):", S$value,sep="")
	}

	newickTree = paste(tree, ';', sep="")
	significativeNodes<-aData$getSignificativeNodes()
	
	if(report){
		report = list()
		report$dims <- c(ncol(data), nrow(data))
		report$initialSimilarityMatrix<-aData$getSimilarityMatrixAtLevel(0)
		report$cor <-cor(data)
		report$freq <- NULL
		report$basicStats <- cbind(apply(data,2,sum), apply(data, 2, mean), apply(data,2, sd))
		colnames(report$basicStats) <- c("freq","mean","sd")
		report$tree <- newickTree
		report$significativeNodes <- significativeNodes
		return(report)
	}
	return(list(tree = newickTree, significativeNodes= significativeNodes))
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