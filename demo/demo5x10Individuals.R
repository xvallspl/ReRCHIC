#' Performs the similarity analysis on a 10 individuals and 5 binary variables dataset.
#' 
#' Simple case with 10 individuals and 5 randomly generated binary variables, not using
#' the algorithm function (yet).
#' 
#' @author Xavier Valls \email{xaviervallspla@@gmail.com}
source("~/master/ReRCHIC/R/ASImodel.R")
data(demo5_10)
aData <- ASImodel$new(data)
tree<-array("",ncol(data)-1)
treeWithDistances<-array("",ncol(data)-1)
for( i in 1:( ncol(data)-1))
{
	S <- aData$getMaximumSimilarity()
	aData$joinClasses(S$nodes)
	genericPair<- aData$findGenericPair(i, S$nodes)
	aData$setGenericImplications(genericPair$pos, i)
	aData$setTypicality(genericPair, at.level=i)
	aData$setContribution(genericPair, at.level=i)
	S$nodes[S$nodes>ncol(data)] = tree[S$nodes[S$nodes>ncol(data)]-ncol(data)]
	print(S$nodes)
	names<-names(S$nodes)
	#tree[i]<-paste("(", S$nodes[1], ",", S$nodes[2],"):", S$value,sep="")
	tree[i]<- paste(colnames(aData$similarityMatrix)[aData$nPrimitiveClasses+i],":", S$value,sep="")
}
treeWithDistances <- colnames(aData$similarityMatrix)[aData$nPrimitiveClasses+1:(ncol(data)-1)]
newickTree = paste(tree, ';', sep="")
newickTreeWithDistances = paste(treeWithDistances, ';', sep="")
print(newickTree)
print(newickTreeWithDistances)
significativeNodes <- aData$getSignificativeNodes()
print(significativeNodes)
