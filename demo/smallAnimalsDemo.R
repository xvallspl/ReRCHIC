#' Performs the similarity analysis on the animalesSmall dataset, plotting the hierarchical tree as a phylogenetic tree.
#' 
#' Similarity analysis on the animalesSmall dataset, a (cropped) dataset describing children perception of animals characteristics.
#' 
#' @author Xavier Valls \email{xaviervallspla@@gmail.com}

data(animalesSmall)
require(ape)
ASI<-callASIAlgorithm(animales)
phylo<-read.tree(text=ASI$tree[length(ASI$tree)])
plot(phylo, use.edge.length=FALSE, direction="upwards")
