report <- setRefClass("reportClass",
	fields = c("data", "initialSimilarityMatrix","bivariantDataParameters", "dataParameters","similarityAtLevel"),
	methods = list( 
			
			initialize = function(ExternalData)
			{	
				# dataParameters <- list(
				# 	frequency = NULL,
				# 	mean = NULL,
				# 	deviation = NULL
				# );

				# bivariantDataParameters <- list( 	
				# 	frequency = NULL,
				# 	mean = NULL,
				# 	deviation = NULL
				# )

				initialSimilarityMatrix <<- NULL
			}
	)
)