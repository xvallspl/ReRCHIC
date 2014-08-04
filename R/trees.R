hierachyTree <- setRefClass("hierarchyTree",
	fields = c("elements", "significative","genericPairs"),
	methods = list( 
			appendHierarchyLevel = function(Tuple,level){
				elements <<- rbind(.self$elements, Tuple)
				Tuple[Tuple>ncol(data)]= tree[Tuple[Tuple>ncol(data)]-ncol(data)]
				tree[level]=paste("(", Tuple[1],",",Tuple[2] ,")", sep="")
				setTipycallity(Tuple)
				setContribution(Tuple)
				calculaIndiceParaSignificativo(Tuple)
			},
			show = function(){},
			getLevel = function(level){},
			setTipycallity = function(Tuple){},
			setContribution = function(Tuple) {},
			calculaIndiceParaSignificativo = function(Tuple){}
	)
)

similarityTree <- setRefClass("similarityTree", contains = "hierarchyTree")

cohesionTree <- setRefClass("cohesionTree", contains = "hierarchyTree")
