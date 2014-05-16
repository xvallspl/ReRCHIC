hierachyTree <- setRefClass("hierarchyTree",
	fields = c("elements", "significative","genericPairs"),
	methods = list( 
			appendHierarchyLevel = function(Tuple){
				elements <<- rbind(.self$elements, Tuple)
				setTipycallity(Tuple)
				setContribution(Tuple)
				calculaIndiceParaSignificativo(Tuple)
			},
			setTipycallity = function(Tuple){},
			setContribution = function(Tuple) {},
			calculaIndiceParaSignificativo = function(Tuple){}
	)
)

similarityTree <- setRefClass("similarityTree", contains = "hierarchyTree")

cohesionTree <- setRefClass("cohesionTree", contains = "hierarchyTree")
