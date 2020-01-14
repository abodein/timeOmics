context("unscale")

test_that("unscale works", {
    
    X <- matrix(1:9, ncol = 3)
    X.scale <- scale(X, center = TRUE, scale = TRUE)
    X.unscale <- unscale(X.scale)
    
    expect_is(X.unscale, "matrix")
    expect_equal(X, X.unscale)
})