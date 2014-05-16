ASIdata <- setRefClass("data",
	fields = c("data", "model", "joinMatrix", "similarityMatrix","genericImplicationsMatrix", "nPrimitiveClasses"),
	methods = list( 
			
			initialize = function( ExternalData, model = "pois" )
			{
				if(!is.matrix(ExternalData))
					 stop("Provided data should be a matrix")
				if(!model %in% c("pois", "binom")) 
					stop("Specified probability distribution model not implemented")
				
				data <<- ExternalData
				nPrimitiveClasses <<- ncol(data)
				.self$model <<- model
				initializeSimilarityMatrix()
				joinMatrix <<- matrix(NaN, nPrimitiveClasses, nPrimitiveClasses)
				diag(joinMatrix)<<-(1:nPrimitiveClasses)
				genericImplicationsMatrix <<- matrix(1, nrow(data), nPrimitiveClasses-1)
			},	

			initializeSimilarityMatrix = function(){
			  	nBinaryPresences <- apply(data,2,'sum')
				nBinaryCopresences <- crossprod(data, data)
			  	probabilityOfCopresence <- tcrossprod(nBinaryPresences, nBinaryPresences)/nPrimitiveClasses^2 

				if(model == 'pois'){
					similarityMatrix <<- ppois( q = nBinaryCopresences, lambda = nPrimitiveClasses*probabilityOfCopresence )
				}else if(model == 'binom'){
					similarityMatrix <<- pbinom( q = nBinaryCopresences, size = nPrimitiveClasses, prob = probabilityOfCopresence )  
				}
				diag(similarityMatrix)<<-0
			},

			getMaximumSimilarity = function(level = ncol(similarityMatrix)- nPrimitiveClasses){
				joined <- getJoinedClassesAtLevel(level)
				if(!length(joined)) 
					simMat <- similarityMatrix
				else 
					simMat <- similarityMatrix[-joined, -joined] 
			 	return(which( simMat == max(simMat), arr.ind = TRUE )[1, ]) 
			 },

			 getJoinedClassesAtLevel = function(level, primitivesOnly = FALSE){
			 	if(level> ncol(similarityMatrix))
					stop("Wrong level or level still not computed")
				if(!primitivesOnly) n <- nPrimitiveClasses+level else n <- nPrimitiveClasses
				return(unique(which(joinMatrix[1:n, 1:n] <= level , arr.ind = T))[, 1])
			 },

			joinClasses = function(Tuple){
			 	computeNewClassSimilarities(Tuple)
			 	updateJoinMatrix(Tuple)
			},

			computeNewClassSimilarities = function(Tuple){
				newClassRow <- newClassCol <- array(0,ncol(similarityMatrix))
				joinedWithFirst  <- which(!(joinMatrix[Tuple[1], ] %in% NaN))
				joinedWithSecond <- which(!(joinMatrix[Tuple[2], ] %in% NaN))
				joinedWithTuple <-c(joinedWithFirst,joinedWithSecond)
				primitivesJoinedWithTuple <- joinedWithTuple[joinedWithTuple<(nPrimitiveClasses+1)]
				for( i in 1:ncol(similarityMatrix)){
					if(!(i %in% joinedWithTuple))
					{	
						classOutsideTuple <- which(!(joinMatrix[i, 1:nPrimitiveClasses] %in% NaN))
						newClassCol[i] <- max(similarityMatrix[classOutsideTuple, primitivesJoinedWithTuple])^length(primitivesJoinedWithTuple)*length(classOutsideTuple)
						newClassRow[i] <- max(similarityMatrix[primitivesJoinedWithTuple, classOutsideTuple])^length(primitivesJoinedWithTuple)*length(classOutsideTuple)
					}
				}
				similarityMatrix <<- rbind(cbind(similarityMatrix, newClassCol), c(newClassRow,0))
			},
			
			updateJoinMatrix = function(Tuple){
				newClass <- array(NaN, ncol(joinMatrix))
				joinedWithFirst	 <- which( !(joinMatrix[Tuple[1], ] %in% NaN))
				joinedWithSecond <- which( !(joinMatrix[Tuple[2], ] %in% NaN))
				joinMatrix[joinedWithFirst,joinedWithSecond] <<- ncol(similarityMatrix)
				joinMatrix[joinedWithSecond,joinedWithFirst] <<- ncol(similarityMatrix)	#Simmetry
				newClass[c(joinedWithFirst, joinedWithSecond)] <- ncol(similarityMatrix)
				joinMatrix <<- rbind(cbind(joinMatrix, newClass), c(newClass, ncol(similarityMatrix)))
			},

			findGenericPairAtLevel = function(Tuple, level=ncol(similarityMatrix)){
				joinedWithFirst  <- which(!(joinMatrix[Tuple[1],] %in% NaN))
				joinedWithSecond <- which(!(joinMatrix[Tuple[2],] %in% NaN))
				
				joinedWithFirst  <- joinedWithFirst[joinedWithFirst < level]
				joinedWithSecond <- joinedWithSecond[joinedWithSecond < level]

				maxPhi <- max(similarityMatrix[joinedWithFirst, joinedWithSecond]) 
				maxPhiInd <- which(similarityMatrix == maxPhi, arr.ind=T)[1,]
				
				return(list(phi = maxPhi, pos = maxPhiInd))
			},

			implicacionesGenericas = function(genericPair, p = 0.5){
				
				for(i in 1:nrow(data))
				{
					if(data[i, genericPair[2]]!=1)
						if(data[i, genericPair[1]]) genericImplicationsMatrix[i, level]<<-0
						else genericImplicationsMatrix[i, level]<<-p
				}
			},

			distance2 = function(individual, class, genericPair){
				return( 1/nrow(data)* sum((genericPair$phi-genericImplicationsMatrix[,class])^2
								  		  / 1-genericPair$phi))
			},

			getSignificativeNodes = function(){
				centeredIndices <- rep(NA, ncol(similarityMatrix)-nPrimitiveClasses)
				for( i in nPrimitiveClasses:ncol(similarityMatrix))
				{	c <- computeCardinalAtLevel(i)
					centeredIndices[i]  <- (c$cardinal- 1/2 * c$nSep * c$nJoin)/
										sqrt(c$nJoin*c$nSep(c$njoin+c$nSep+1)/12)
				}
				v <- diff(diff(centeredIndex))
				localMaximumPositions <- c(1, which(diff(sign(v))!=0))
				return(localMaximumPositions)
			},

			computeCardinalAtLevel = function(level){
				primJoined <- getJoinedClassesAtLevel(level, primitivesOnly = TRUE)
				primSepparated <- (1:nPrimitiveClasses)[-primaryJoined]
				
				SimilaritiesJoined     <- similarityMatrix[primJoined, primJoined]
				SimilaritiesSepparated <- similarityMatrix[primSepparated, primSepparated]
				# Mann - Whitney
				cardinal <- sum(which(sort(unique(c(SimilaritiesJoined, SimilaritiesSepparated))) 
										  %in% SimilaritiesSepparated)) - 
				            (length(SimilaritiesJoined)*(length(SimilaritiesJoined)+1))/(2*2) 
				return(list( cardinal= cardinal, nJoin = length(SimilaritiesJoined)/2, nSep = length(SimilaritiesSepparated)/2))
			}
	)
)
