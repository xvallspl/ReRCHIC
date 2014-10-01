#' A Reference Class to represent the Statistical Implicative Analysis model.
#' @name ASImodel
# @import methods
#' @export ASImodel
#' @exportClass ASImodel
#' @field data   External data.
#' @field model Probability model to use. Binom or pois.
#' @field joinMatrix Matrix of joint classes
#' @field similarityMatrix Matrix of bivariant similarities.
#' @field genericImplicationsMatrix Matrix of generic Implications
#' @field nPrimitiveClasses Number of primitive classes

ASImodel <- setRefClass("ASImodel",
  fields = c("data", "model", "levels", "joinMatrix", "similarityMatrix",
  			 "genericImplicationsMatrix", "joinedClasses","nPrimitiveClasses", "contributionMatrix",
  			 "typicalityMatrix","nIndividuals"),
  methods = list( 

      initialize = function( externalData, model = "pois" )
      {
        if(!is.matrix(externalData))
           stop("Provided data should be a matrix")
        if(!model %in% c("pois", "binom")) 
          stop("Specified probability distribution model not implemented")
        
        nIndividuals <<- nrow(externalData)
        nPrimitiveClasses <<- ncol(externalData)
        data<<- externalData
        .self$model <<- model
        initializeSimilarityMatrix(externalData)
        joinMatrix <<- matrix(NaN, nPrimitiveClasses, nPrimitiveClasses)
        diag(joinMatrix)<<-(Inf)
        genericImplicationsMatrix <<- matrix(1, nIndividuals, nPrimitiveClasses-1)
        contributionMatrix <<- matrix(NaN, nIndividuals, nPrimitiveClasses-1)
        typicalityMatrix <<- matrix(NaN, nIndividuals, nPrimitiveClasses-1)
        joinedClasses<<-list()
        levels <<- 0
      },  

      initializeSimilarityMatrix = function(externalData){
          nBinaryPresences <- apply(externalData,2,'sum')
        nBinaryCopresences <- crossprod(externalData, externalData)
          probabilityOfCopresence <- tcrossprod(nBinaryPresences, nBinaryPresences)/nIndividuals^2 

        if(model == 'pois'){
          similarityMatrix <<- ppois( q = nBinaryCopresences, 
          							  lambda = nIndividuals*probabilityOfCopresence )
        }else if(model == 'binom'){
          similarityMatrix <<- pbinom( q = nBinaryCopresences, size = nIndividuals, 
          							   prob = probabilityOfCopresence )  
        }
        diag(similarityMatrix)<<-0
      },

      getMaximumSimilarity = function(at.level = levels){
        "Gets the maximum similartity at the specified level. By default the last one"
        simMat <- getSimilarityMatrix(at.level)
        joined <- getJoinedClasses(at.level)
        M <- which( similarityMatrix == max(simMat), arr.ind = TRUE )
        MnotJoined <- apply( M, 1, function(x){all(!(x %in%joined))})
        nodes <- M[MnotJoined,][1, ]
        names(nodes)<-colnames(similarityMatrix)[nodes]
         return(list(nodes = nodes, value = max(simMat)))
       },

       getJoinedClasses = function(at.level=levels, primitivesOnly = FALSE){
         .checkIfLevelExists(at.level)
         if(at.level==0) 
           ret=c()
         else{
           ret<-joinedClasses[[at.level]]
          if(primitivesOnly) ret<-ret[ret<=nPrimitiveClasses]
        }
        return(ret)
       },

       joinedWithClass = function(class, at.level=levels, primitivesOnly = FALSE){
         .checkIfLevelExists(at.level)
        if(!primitivesOnly){
          ret <- class
          n   <- nPrimitiveClasses+at.level 
        }
        else{
          ret <- c()
          n   <- nPrimitiveClasses
        }   
        joined <- which(!(joinMatrix[class, 1:n]  %in% NaN) & 
                     joinMatrix[class, 1:n] <=at.level)
         return(c(ret, joined))
       },

      joinClasses = function(Tuple){
        "Joins two classes into a new node of the hierarchical tree"
        if(ncol(similarityMatrix) == (2*nPrimitiveClasses-1)){ 
        	stop("You're already at the last level!")
        }
         .computeNewClassSimilarities(Tuple)
         .updateJoinMatrix(Tuple)
         levels <<- levels+1
         if(levels>1)
         {
           joinedClasses[[levels]]<<-c(joinedClasses[[levels-1]],Tuple)
         }
         else joinedClasses[[levels]]<<-Tuple
      },

      .computeNewClassSimilarities = function(Tuple){
        newClassRow <- newClassCol <- array(0,ncol(similarityMatrix))
        joinedWithFirst  <- joinedWithClass(Tuple[1])
        joinedWithSecond <- joinedWithClass(Tuple[2])
        joinedWithTuple <-c(joinedWithFirst,joinedWithSecond)
        primitivesJoinedWithTuple <- joinedWithTuple[joinedWithTuple<(nPrimitiveClasses+1)]
        for( i in 1:ncol(similarityMatrix)){
          if(!(i %in% joinedWithTuple))
          {  
            classOutsideTuple <- which(!(joinMatrix[i, 1:ncol(similarityMatrix)] %in% NaN))
            newClassCol[i] <- max(similarityMatrix[classOutsideTuple, primitivesJoinedWithTuple])^
            						(length(primitivesJoinedWithTuple)*length(classOutsideTuple))
            newClassRow[i] <- max(similarityMatrix[primitivesJoinedWithTuple, classOutsideTuple])^
            						(length(primitivesJoinedWithTuple)*length(classOutsideTuple))
          }
        }
        similarityMatrix <<- rbind(cbind(similarityMatrix, newClassCol), c(newClassRow,0))
        colnames(similarityMatrix)[ncol(similarityMatrix)]<<- paste("(",
        		colnames(similarityMatrix)[Tuple[1]],",", colnames(similarityMatrix)[Tuple[2]],")")        
      },
      
      .updateJoinMatrix = function(Tuple){
        newClass <- array(NaN, ncol(joinMatrix))
        joinedWithFirst  <- joinedWithClass(Tuple[1])
        joinedWithSecond <- joinedWithClass(Tuple[2])
        joinMatrix[joinedWithFirst,joinedWithSecond] <<- levels+1
        joinMatrix[joinedWithSecond,joinedWithFirst] <<- levels+1  #Simmetry
        newClass[c(joinedWithFirst, joinedWithSecond)] <- levels+1
        joinMatrix <<- rbind(cbind(joinMatrix, newClass), c(newClass, Inf))
        colnames(joinMatrix)[ncol(joinMatrix)] <<-levels+1
      },

      findGenericPair = function(level, Tuple){
        "Finds the generic Pair of a derived class"
        joinedWithFirst  <- joinedWithClass(Tuple[1], level, primitivesOnly=TRUE)
        joinedWithSecond <- joinedWithClass(Tuple[2], level, primitivesOnly=TRUE)

        maxPhi <- max(similarityMatrix[joinedWithFirst, joinedWithSecond])
        simMat <- similarityMatrix[1:nPrimitiveClasses, 1:nPrimitiveClasses]
        maxPhiInd <- which( simMat == maxPhi, arr.ind=T)[1,]
        
        return(list(phi = maxPhi, pos = maxPhiInd))
      },

      setGenericImplications = function(genericPair, level, p = 0.5){
        "Sets the generic implication of each individual for a given derived class (level)"        
        for(i in 1:nIndividuals)
        {
          if(data[i, genericPair[2]]!=1)
            if(data[i, genericPair[1]]) genericImplicationsMatrix[i, level]<<-0
            else genericImplicationsMatrix[i, level]<<-p
        }
      },

      getSignificativeNodes = function(){
        "Gets the significative nodes of the hierarchical tree"
        centeredIndices <- rep(NA, ncol(similarityMatrix)-nPrimitiveClasses)
        for( i in 1:(ncol(similarityMatrix)-nPrimitiveClasses))
        {  c <- computeCardinal(i)
          centeredIndices[i]  <- (c$cardinal- 1/2 * c$nSep * c$nJoin)/
                    sqrt(c$nJoin*c$nSep*(c$nJoin+c$nSep+1)/12)
        }
        v <- diff(diff(centeredIndices))
        localMaximumPositions <- which(diff(sign(v))!=0)
        return(localMaximumPositions)
      },

      computeCardinal = function(at.level=levels){
        primJoined <- getJoinedClasses(at.level, primitivesOnly = TRUE)
        primSepparated <- (1:nPrimitiveClasses)[-primJoined]
        SimilaritiesJoined     <- similarityMatrix[primJoined, primJoined]
        SimilaritiesSepparated <- similarityMatrix[primSepparated, primSepparated]
        # Mann - Whitney
        cardinal <- sum(which(sort(unique(c(SimilaritiesJoined, SimilaritiesSepparated))) 
                      %in% SimilaritiesSepparated)) - 
                    (length(SimilaritiesJoined)*(length(SimilaritiesJoined)+1))/(2*2) 
        return(list( cardinal= cardinal, nJoin = length(SimilaritiesJoined)/2,
        			 nSep = length(SimilaritiesSepparated)/2))
      },
      
      getSimilarityMatrix = function(at.level=levels){
        "Returns the similarity matrix at the specified level of the hierarchical tree"
        .checkIfLevelExists(at.level)
        omit <-getJoinedClasses(at.level)
        if(!length(omit)) 
          simMat <- similarityMatrix[1:nPrimitiveClasses, 1:nPrimitiveClasses]
        else 
          simMat <- similarityMatrix[1:(nPrimitiveClasses+at.level),
          							 1:(nPrimitiveClasses+at.level)][-omit, -omit] 
        return(simMat)
      },

      .checkIfLevelExists = function(level){
        if((level+nPrimitiveClasses)> ncol(similarityMatrix))
          stop("Wrong level or level still not computed")
      },

      distance2 = function(individual, level, genericPair){
        classSubClasses <- .getSubLevels(level)
        return( 1/length(classSubClasses) * sum((genericPair$phi-
        	genericImplicationsMatrix[individual,classSubClasses])^2 / (1-genericPair$phi)))
      },

      .getSubLevels = function(class){
        joinedWithClass <- joinedWithClass(nPrimitiveClasses+class, class)
        classSubClasses <- joinedWithClass[which(
        		joinMatrix[nPrimitiveClasses+class,joinedWithClass] == class)]
        return(c(class, classSubClasses[classSubClasses > nPrimitiveClasses]-nPrimitiveClasses))
      },

      distanceTilde2 = function(individual, level){
        classSubclasses <- .getSubLevels(level)
        return( 1/length(classSubclasses) * sum((1-genericImplicationsMatrix[individual,
        																	classSubclasses])^2))
      },      

      setTypicality = function(genericPair, at.level=levels ){        
        "Computes the typicality of each individual for a certain class"
        dist<-array(NaN, nIndividuals)
        for(i in 1:nIndividuals)
        {
          dist[i]<-sqrt(distance2(i, at.level, genericPair))
        }
        typicalityMatrix[,at.level]<<- 1 - dist/ max(dist)
      },

      setContribution = function(genericPair, at.level=levels){
        "Computes the contribution of each individual for a certain class"        
        for(i in 1:nIndividuals)
        {
          contributionMatrix[i,at.level]<<-(1 - sqrt(distanceTilde2(i, at.level)))
        }
      }
  )
)
