#' Performs the similarity analysis on arules' Adult dataset, providing a report.
#' 
#' Similarity analysis with report on arules' Adult dataset with a boolean matrix instead of an integer binary one.
#' 
#' @author Xavier Valls \email{xaviervallspla@@gmail.com}

require(arules)
data(Adult)
data<-as.matrix(Adult@data)[,1:100]
report<-callASIAlgorithm(data, report = TRUE)
displayReport(report)
