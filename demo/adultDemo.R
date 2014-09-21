#' Performs the similarity analysis on the animalesSmall dataset, plotting the hierarchical tree as a phylogenetic tree.
#' 
#' Similarity analysis on arules' dataset, a dataset describing children perception of animals characteristics.
#' 
#' @author Xavier Valls \email{xaviervallspla@@gmail.com}

require(arules)
data(Adult)
data<-as.matrix(Adult@data)[,1:100]
report<-callASIAlgorithm(data, report = TRUE)
displayReport(report)
