context("tuneCluster.spca")

demo <- suppressMessages(get_demo_cluster())
X <- demo$X
tune.spca.res <- tuneCluster.spca(X = X, ncomp = 2, test.keepX = c(2,5,7))
# plot(tune.spca.res)
# plot(tune.spca.res, comp = 2)

test_that("tuneCluster.spca failed on invalid input - X", {
    expect_error(tuneCluster.spca(X = ""), "X must be a numeric matrix/data.frame", fixed = TRUE)
    expect_error(tuneCluster.spca(X = 1), "X must be a numeric matrix/data.frame", fixed = TRUE)
    expect_error(tuneCluster.spca(X = NA), "X must be a numeric matrix/data.frame", fixed = TRUE)
    expect_error(tuneCluster.spca(X = list()), "X must be a numeric matrix/data.frame", fixed = TRUE)
})

test_that("tuneCluster.spca failed on invalid input - ncomp", {
    expect_error(tuneCluster.spca(X = demo$X, ncomp = ""), "'ncomp' should be an integer between 1 and 10", fixed = TRUE)
    expect_error(tuneCluster.spca(X = demo$X, ncomp = c(1,2)), "'ncomp' should be an integer between 1 and 10", fixed = TRUE)
    expect_error(tuneCluster.spca(X = demo$X, ncomp = 0), "'ncomp' should be an integer between 1 and 10", fixed = TRUE)
    # ncomp < nrow(X)
    expect_error(tuneCluster.spca(X = demo$X, ncomp = 11),"'ncomp' should be an integer between 1 and 10", fixed = TRUE)
})

test_that("tuneCluster.spca failed on invalid input - keepX", {
    expect_error(tuneCluster.spca(X = demo$X, ncomp = 2, test.keepX = c("a",1)), "'test.keepX' should be numeric", fixed = TRUE)
    expect_error(tuneCluster.spca(X = demo$X, ncomp = 2, test.keepX = list()), "'test.keepX' should be numeric", fixed = TRUE)
    expect_error(tuneCluster.spca(X = demo$X, ncomp = 2, test.keepX = "abc"), "'test.keepX' should be numeric", fixed = TRUE)
    expect_error(tuneCluster.spca(X = demo$X, ncomp = 2, test.keepX = matrix(1:9)), "'test.keepX' should be numeric", fixed = TRUE)
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
