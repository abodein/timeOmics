context("tuneCluster.spca")

demo <- suppressMessages(get_demo_cluster())
X <- demo$X
tune.spca.res <- tuneCluster.spca(X = X, ncomp = 2, test.keepX = c(2,5,7))
# plot(tune.spca.res)
# plot(tune.spca.res, comp = 2)

test_that("tuneCluster.spca failed on invalid input", {
    #-- X should be a numeric matrix/data.frame
    expect_error(tuneCluster.spca(X = ""), "X must be a numeric matrix with finite value", fixed = TRUE)
    expect_error(tuneCluster.spca(X = 1), "X must be a numeric matrix with finite value", fixed = TRUE)
    expect_error(tuneCluster.spca(X = NA), "X must be a numeric matrix with finite value", fixed = TRUE)
    expect_error(tuneCluster.spca(X = list()), "X must be a numeric matrix with finite value", fixed = TRUE)
    
    #-- ncomp
    expect_error(tuneCluster.spca(X = demo$X, ncomp = ""), "invalid value for 'ncomp'", fixed = TRUE)
    expect_error(tuneCluster.spca(X = demo$X, ncomp = c(1,2)), "invalid value for 'ncomp'", fixed = TRUE)
    expect_error(tuneCluster.spca(X = demo$X, ncomp = 0), "invalid value for 'ncomp'", fixed = TRUE)
    # ncomp < nrow(X)
    expect_error(tuneCluster.spca(X = demo$X, ncomp = 11), "use smaller 'ncomp'", fixed = TRUE)
    
    #-- keepX
    expect_error(tuneCluster.spca(X = demo$X, ncomp = 2, test.keepX = c("a",1)), "'test.keepX' must be a numeric vector with more than two entries", fixed = TRUE)
    expect_error(tuneCluster.spca(X = demo$X, ncomp = 2, test.keepX = 1), "'test.keepX' must be a numeric vector with more than two entries", fixed = TRUE)
    expect_error(tuneCluster.spca(X = demo$X, ncomp = 2, test.keepX = list()), "'test.keepX' must be a numeric vector with more than two entries", fixed = TRUE)
    expect_error(tuneCluster.spca(X = demo$X, ncomp = 2, test.keepX = "abc"), "'test.keepX' must be a numeric vector with more than two entries", fixed = TRUE)
    expect_error(tuneCluster.spca(X = demo$X, ncomp = 2, test.keepX = matrix(1:9)), "'test.keepX' must be a numeric vector with more than two entries", fixed = TRUE)
    expect_error(tuneCluster.spca(X = demo$X, ncomp = 2, test.keepX = NULL), "'test.keepX' must be a numeric vector with more than two entries", fixed = TRUE)
})

test_that("tuneCluster.spca works", {
    expect_is(tuneCluster.spca(X = X, ncomp = 2, test.keepX = c(5:10)), "spca.tune.silhouette")
    expect_is(tuneCluster.spca(X = as.data.frame(X), ncomp = 3, test.keepX = c(5:10)), "spca.tune.silhouette")
})

test_that("plot.spca.tune.silhouette failed on invalid input", {
    #-- comp
    expect_error(plot(tune.spca.res, comp = 3), "Invalid 'comp', shoud be an integer between 1 and 2")
    expect_error(plot(tune.spca.res, comp = NA),"Invalid 'comp', shoud be an integer between 1 and 2")
    expect_error(plot(tune.spca.res, comp = c(1,2)),"Invalid 'comp', shoud be an integer between 1 and 2")
})

test_that("plot.spca.tune.silhouette works", {
    #-- comp
    expect_is(plot(tune.spca.res, comp = 1), "ggplot")
    expect_is(plot(tune.spca.res, comp = 2), "ggplot")
    expect_is(plot(tune.spca.res, comp = 2, plot = TRUE), "ggplot")
    expect_is(plot(tune.spca.res, comp = 2, plot = FALSE), "ggplot")
    expect_is(plot(tune.spca.res, comp = 2, plot = NULL), "ggplot")
})
