context("tuneCluster.spls")

demo <- get_demo_cluster()
X <- demo$X
Y <- demo$Y
res <- spls(X, Y, ncomp = 10)
#tune.spls <- tuneCluster.spls(X, Y, ncomp = 2, test.keepX = c(5,10,15,20), test.keepY <- c(2,4,6))

test_that("tuneCluster.spls failed on invalid input - X", {
    #-- X must be a numeric matrix/data.frame
    expect_error(tuneCluster.spls(X = ""), "X must be a numeric matrix/data.frame")
    expect_error(tuneCluster.spls(X = list()), "X must be a numeric matrix/data.frame")
    expect_error(tuneCluster.spls(X = NULL), "X must be a numeric matrix/data.frame")
    expect_error(tuneCluster.spls(X = NA), "X must be a numeric matrix/data.frame")
    expect_error(tuneCluster.spls(X = matrix(letters[1:9], ncol = 3)), "X must be a numeric matrix/data.frame")
})

test_that("tuneCluster.spls failed on invalid input - X", {
    #-- X must be a numeric matrix/data.frame
    expect_error(tuneCluster.spls(X = demo$X, Y = ""), "Y must be a numeric matrix/data.frame")
    expect_error(tuneCluster.spls(X = demo$X, Y = list()), "Y must be a numeric matrix/data.frame")
    expect_error(tuneCluster.spls(X = demo$X, Y = NULL), "Y must be a numeric matrix/data.frame")
    expect_error(tuneCluster.spls(X = demo$X, Y = NA), "Y must be a numeric matrix/data.frame")
    expect_error(tuneCluster.spls(X = demo$X, Y = matrix(letters[1:9], ncol = 3)), "Y must be a numeric matrix/data.frame")
})


test_that("tuneCluster.spls failed on invalid input - ncomp", {   
    #-- ncomp  // integer 
    expect_error(tuneCluster.spls(X = X, Y = Y, ncomp = NULL), "'ncomp' should be an integer between 1 and 10", fixed = TRUE)
    expect_error(tuneCluster.spls(X = X, Y = Y, ncomp = "abc"), "'ncomp' should be an integer between 1 and 10", fixed = TRUE)
    expect_error(tuneCluster.spls(X = X, Y = Y, ncomp = X), "'ncomp' should be an integer between 1 and 10", fixed = TRUE)
    expect_error(tuneCluster.spls(X = X, Y = Y, ncomp = 55), "'ncomp' should be an integer between 1 and 10", fixed = TRUE)
})

test_that("tuneCluster.spls failed on invalid input - test.keepX", {  
    #-- test.keepX  // numeric vecter lower than ncol(X)
    expect_error(tuneCluster.spls(X = X, Y = Y, test.keepX = "abc"), "'test.keepX' should be numeric")
    expect_error(tuneCluster.spls(X = X, Y = Y, test.keepX = NA), "'test.keepX' should be numeric")
    expect_error(tuneCluster.spls(X = X, Y = Y, test.keepX = matrix()), "'test.keepX' should be numeric")
    expect_error(tuneCluster.spls(X = X, Y = Y, test.keepX = list()), "'test.keepX' should be numeric")
    expect_error(tuneCluster.spls(X = X, Y = Y, test.keepX = c(1,2,101)), "'test.keepX' must be lower than 100, ncol(X)", fixed=TRUE)
})

test_that("tuneCluster.spls failed on invalid input - test.keepY", {  
    #-- test.keepX  // numeric vecter lower than ncol(X)
    expect_error(tuneCluster.spls(X = X, Y = Y, test.keepY = "abc"), "'test.keepY' should be numeric")
    expect_error(tuneCluster.spls(X = X, Y = Y, test.keepY = NA), "'test.keepY' should be numeric")
    expect_error(tuneCluster.spls(X = X, Y = Y, test.keepY = matrix()), "'test.keepY' should be numeric")
    expect_error(tuneCluster.spls(X = X, Y = Y, test.keepY = list()), "'test.keepY' should be numeric")
    expect_error(tuneCluster.spls(X = X, Y = Y, test.keepY = c(1,2,101)), "'test.keepY' must be lower than 10, ncol(Y)", fixed=TRUE)
})

test_that("tuneCluster.spls works", {  
    #-- test.keepX  // numeric vecter lower than ncol(X)
    expect_is(tuneCluster.spls(X = X, Y = Y), "spls.tune.silhouette")
    expect_is(tuneCluster.spls(X = X, Y = Y, ncomp = 3), "spls.tune.silhouette")
    expect_is(tuneCluster.spls(X = X, Y = Y, ncomp = 3, test.keepX = NULL, test.keepY = NULL), "spls.tune.silhouette")
    expect_is(tuneCluster.spls(X = X, Y = Y, ncomp = 2, test.keepX = c(1,2)), "spls.tune.silhouette")
})