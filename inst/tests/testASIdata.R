context("ASImodel")

externalData <- matrix(0, 5, 5)
data <- as.matrix(read.csv2(file="test-example.csv", header=TRUE, row.names=1))
aData <- ASImodel$new(data)

test_that("If provided data isn't a matrix raises error",{
	expect_error(ASImodel$new(c(1,2)), "Provided data should be a matrix")
	expect_error(ASImodel$new(1) , "Provided data should be a matrix")
})

test_that("checks provided models",
{
	expect_error(ASImodel$new(externalData, "student t"), "Specified probability distribution model not implemented")
	expect_error(ASImodel$new(externalData, ""), "Specified probability distribution model not implemented")
})

test_that("ASImodel initializes on creation",{
	aData <- ASImodel$new(externalData)
	expect_that(aData$data, equals(externalData))
	expect_that(aData$nPrimitiveClasses, equals(5))
	expect_that(aData$model, equals("pois"))
})

test_that("Obtains correctly the joined classes at the specified level",{
	Tuple <- aData$getMaximumSimilarity()
	aData$joinClasses(Tuple)
	joined <- aData$getJoinedClassesAtLevel(1)
	expect_equal(c(Tuple,ncol(aData$getSimilarityMatrixAtLevel(1))), joined)
})

# test_that("",{

# })

# 			computeSimilarityMatrix = function(){
			
# 			getMaximumSimilarity = function(level = ncol(similarityMatrix)){
			
# 			getJoinedClassesAtLevel = function(level, prim = FALSE){

# 			joinClasses = function(Tuple){

# 			buildNewClass = function(Tuple){

# 			markJoinsInJoinMatrix = function(Tuple){

# 			getSignificativeNodes = function(){

# 			computeCardinalAtLevel = function(level){

# 			findGenericPairAtLevel = function(Tuple, level=ncol(data)){


# test_that("buildNewClass returns the correct class",{
# 	expect_that( aData$buildNewClass, )
# })

