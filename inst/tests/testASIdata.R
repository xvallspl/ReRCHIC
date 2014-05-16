context("ASIdata")

externalData <- matrix(0, 5, 5)

test_that("If provided data isn't a matrix raises error",{
	expect_error(ASIdata$new(c(1,2)), "Provided data should be a matrix")
	expect_error(ASIdata$new(1) , "Provided data should be a matrix")
})

test_that("checks provided models",
{
	expect_error(ASIdata$new(externalData, "student t"), "Specified probability distribution model not implemented")
	expect_error(ASIdata$new(externalData, ""), "Specified probability distribution model not implemented")
})

test_that("ASIdata initializes on creation",{
	aData <- ASIdata$new(externalData)
	expect_that(aData$data, equals(externalData))
	expect_that(aData$indices, equals(1:5))
	expect_that(aData$nPrimaryClasses, equals(5))
	expect_that(aData$model, equals("pois"))
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

