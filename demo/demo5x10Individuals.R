#' Performs the similarity analysis on a 10 individuals and 5 binary variables dataset.
#' 
#' Simple case with 10 individuals and 5 randomly generated binary variables, not using
#' the algorithm function (yet).
#' 
#' @author Xavier Valls \email{xaviervallspla@@gmail.com}

data(demo5_10)
aData <- ASImodel$new(data)
tree<-array("",ncol(data)-1)
for( i in 1:( ncol(data)-1))
	{
		S <- aData$getMaximumSimilarity()
		#append(tree$genericPairs, aData$findGenericPairAtLevel(S$nodes))
		aData$joinClasses(S$nodes)
		S$nodes[S$nodes>ncol(data)] = tree[S$nodes[S$nodes>ncol(data)]-ncol(data)]
		tree[i]=paste("(", S$nodes[1],",",S$nodes[2] ,"):", S$value,sep="")
	}
	newickTree = paste(tree, ';', sep="")
report = list()
report$initialSimilarityMatrix<-aData$getSimilarityMatrixAtLevel(0)
report$cor <-cor(data)
report$basicStats <- cbind(apply(data,2,sum), apply(data, 2, mean), apply(data,2, sd))
colnames(report$basicStats)<-c("freq","mean","sd")
