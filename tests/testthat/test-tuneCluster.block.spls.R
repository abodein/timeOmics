context("tuneCluster.block.spls")

demo <- suppressMessages(get_demo_cluster())
X <- list(X = demo$X, Z = demo$Z)
Y <- demo$Y
test.list.keepX <- list("X" = c(5,10,15), "Z" = c(2,4,6))
test.keepY <- c(2:5)

#tune.block.spls <- tuneCluster.block.spls(X = X, Y = Y, test.list.keepX = test.list.keepX, test.keepY = test.keepY, mode = "canonical")

test_that("tuneCluster.block.spls failed on invalid input - X", {
    #-- X  // list of data.frame/matrix
    expect_error(tuneCluster.block.spls(X = demo$X), "X must be a list of matrix/data.frame", fixed = TRUE)
    expect_error(tuneCluster.block.spls(X = NA), "X must be a list of matrix/data.frame", fixed = TRUE)
    expect_error(tuneCluster.block.spls(X = NULL), "X must be a list of matrix/data.frame", fixed = TRUE)
    expect_error(tuneCluster.block.spls(X = "ah!"), "X must be a list of matrix/data.frame", fixed = TRUE)
}) 

test_that("tuneCluster.block.spls failed on invalid input - Y", {   
    #-- Y  // NULL or matrix
    expect_error(tuneCluster.block.spls(X = X, Y = ""), "Y must be a numeric matrix/data.frame", fixed = TRUE)
    expect_error(tuneCluster.block.spls(X = X, Y = NA), "Y must be a numeric matrix/data.frame", fixed = TRUE)
}) 

test_that("tuneCluster.block.spls failed on invalid input - indY", {   
    #-- indY  // if Y is NULL, numeric
    expect_error(tuneCluster.block.spls(X = X, Y = NULL), "'indY' must be a numeric value lower or equal to 2, the number of blocks in X.", fixed = TRUE)
    expect_error(tuneCluster.block.spls(X = X, indY = 3), "'indY' must be a numeric value lower or equal to 2, the number of blocks in X.", fixed = TRUE)
    expect_error(tuneCluster.block.spls(X = X, indY = ""), "'indY' must be a numeric value lower or equal to 2, the number of blocks in X.", fixed = TRUE)
    expect_error(tuneCluster.block.spls(X = X, indY = NULL), "'indY' must be a numeric value lower or equal to 2, the number of blocks in X.", fixed = TRUE)
}) 

test_that("tuneCluster.block.spls failed on invalid input - ncomp", {   
    #-- ncomp  // integer 
    expect_error(tuneCluster.block.spls(X = X, Y = Y, ncomp = NULL), "'ncomp' should be an integer between 1 and 10", fixed = TRUE)
    expect_error(tuneCluster.block.spls(X = X, Y = Y, ncomp = "abc"), "'ncomp' should be an integer between 1 and 10", fixed = TRUE)
    expect_error(tuneCluster.block.spls(X = X, Y = Y, ncomp = X), "'ncomp' should be an integer between 1 and 10", fixed = TRUE)
    expect_error(tuneCluster.block.spls(X = X, Y = Y, ncomp = 55), "'ncomp' should be an integer between 1 and 10", fixed = TRUE)
}) 

test_that("tuneCluster.block.spls failed on invalid input - test.list.keepX", {     
    #-- test.list.keepX  // list of integer of the same size as X
    expect_error(tuneCluster.block.spls(X = X, Y = Y), "'test.list.keepX' must be a list of numeric of size 2.", fixed = T)
    expect_error(tuneCluster.block.spls(X = X, Y = Y, test.list.keepX = c()), "'test.list.keepX' must be a list of numeric of size 2.", fixed = T)
    expect_error(tuneCluster.block.spls(X = X, Y = Y, test.list.keepX = list()), "'test.list.keepX' must be a list of numeric of size 2.", fixed = T)
    expect_error(tuneCluster.block.spls(X = X, Y = Y, test.list.keepX = list("A" = 1:3, "B"= c())), "'list.test.keepX' should have the same names as X", fixed = T)
    expect_error(tuneCluster.block.spls(X = X, Y = Y, test.list.keepX = list("X" = 1:3, "Z"= c(4,100))), "'test.list.keepX[[2]] sould be lower than 50, ncol(X[[2]]).", fixed = T)
}) 

test_that("tuneCluster.block.spls failed on invalid input - test.keepY", {      
    #-- test.keepY  // vector of size 
    expect_error(tuneCluster.block.spls(X = X, Y = Y, test.list.keepX = list("X" = 1:3, "Z"= c(4,5)), test.keepY = ""), "'test.keepY' should be numeric")
    expect_error(tuneCluster.block.spls(X = X, Y = Y, test.list.keepX = list("X" = 1:3, "Z"= c(4,5)), test.keepY = list()), "'test.keepY' should be numeric")
    expect_error(tuneCluster.block.spls(X = X, Y = Y, test.list.keepX = list("X" = 1:3, "Z"= c(4,5)), test.keepY = NA), "'test.keepY' should be numeric")
})
    
test_that("tuneCluster.block.spls works", {
    expect_is(suppressWarnings(tuneCluster.block.spls(X=X, Y=Y, test.list.keepX = list("X" = 1:3, "Z"= c(4,5)), test.keepY = c(5))),"block.spls.tune.silhouette")
    expect_is(suppressWarnings(tuneCluster.block.spls(X=X, Y=Y, ncomp = 3, test.list.keepX = list("X" = 3, "Z"= c(4,5)), test.keepY = c(5))),"block.spls.tune.silhouette")
    expect_is(suppressWarnings(tuneCluster.block.spls(X=list(X=demo$X, Y=demo$Y, Z=demo$Z), indY=2, ncomp = 2, test.list.keepX = list("X" = 1:3, "Z"= c(4,5), "Y"= 2))),"block.spls.tune.silhouette")
})
