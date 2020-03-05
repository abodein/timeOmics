context("getNcomp")

demo <- get_demo_cluster()
demo$pca <- mixOmics::pca(X = demo$X, ncomp = 5)
demo$pls <- mixOmics::pls(X = demo$X, Y=demo$Y, ncomp =5, mode = "canonical")
demo$block.pls = mixOmics::block.pls(X=list(X=demo$X, Y=demo$Y, Z=demo$Z), indY=1, ncomp = 5, mode = "canonical")


test_that("getNcomp failed on invalid input - object", {
    # test for "object"
    lapply(list("",1,demo$X, NA), function(i) expect_error(getNcomp(i), "invalid object"))
})

test_that("getNcomp failed on invalid input - max.ncomp", {
    expect_error(getNcomp(object=demo$pca, max.ncomp = 0), "'max.ncomp' should be greater than 1")
    expect_error(getNcomp(object=demo$block.pls, max.ncomp = 0), "'max.ncomp' should be greater than 1")
    expect_error(getNcomp(object=demo$block.pls, max.ncomp = 11), "use smaller 'max.ncomp'")
    expect_error(getNcomp(object=demo$pca, max.ncomp = 11), "use smaller 'max.ncomp'")
})

test_that("getNcomp failed on invalid input - X", {
    # pca
    expect_error(getNcomp(object=demo$pca, X = NA), "X must be a numeric matrix/data.frame", fixed = TRUE)
    expect_error(getNcomp(object=demo$pca, X = NULL), "X must be a numeric matrix/data.frame", fixed = TRUE)
    expect_error(getNcomp(object=demo$pca, X = "ah!"), "X must be a numeric matrix/data.frame", fixed = TRUE)
    
    # block.pls
    expect_error(getNcomp(object=demo$block.pls, X = demo$X), "X must be a list of matrix/data.frame", fixed = TRUE)
    expect_error(getNcomp(object=demo$block.pls, X = NA), "X must be a list of matrix/data.frame", fixed = TRUE)
    expect_error(getNcomp(object=demo$block.pls, X = NULL), "X must be a list of matrix/data.frame", fixed = TRUE)
    expect_error(getNcomp(object=demo$block.pls, X = "ah!"), "X must be a list of matrix/data.frame", fixed = TRUE)
})
     
test_that("getNcomp failed on invalid input - Y", {
    #-- Y  // NULL or matrix
    expect_error(getNcomp(object=demo$block.pls, X = list(X =demo$X, Z=demo$Z), Y = ""), "Y must be a numeric matrix/data.frame", fixed = TRUE)
    expect_error(getNcomp(object=demo$block.pls, X = list(X =demo$X, Z=demo$Z), Y = NA), "Y must be a numeric matrix/data.frame", fixed = TRUE)
})

test_that("getNcomp failed on invalid input - indY", {   
    #-- indY  // if Y is NULL, numeric
    expect_error(getNcomp(object=demo$block.pls, X = list(X =demo$X, Z=demo$Z), Y = NULL), "'indY' must be a numeric value lower or equal to 2, the number of blocks in X.", fixed = TRUE)
    expect_error(getNcomp(object=demo$block.pls, X = list(X =demo$X, Z=demo$Z), indY = 3), "'indY' must be a numeric value lower or equal to 2, the number of blocks in X.", fixed = TRUE)
    expect_error(getNcomp(object=demo$block.pls, X = list(X =demo$X, Z=demo$Z), indY = ""), "'indY' must be a numeric value lower or equal to 2, the number of blocks in X.", fixed = TRUE)
    expect_error(getNcomp(object=demo$block.pls, X = list(X =demo$X, Z=demo$Z), indY = NULL), "'indY' must be a numeric value lower or equal to 2, the number of blocks in X.", fixed = TRUE)
}) 

test_that("getNcomp failed works", {  
    expect_is(getNcomp(demo$pca, max.ncomp = 4, X = demo$X), "ncomp.tune.silhouette")
    expect_is(getNcomp(demo$pca, max.ncomp = 4, X = demo$X, scale = TRUE, center = TRUE), "ncomp.tune.silhouette")
    expect_is(getNcomp(demo$pls, max.ncomp = 4, X = demo$X, Y=demo$Y, scale = TRUE), "ncomp.tune.silhouette")
    expect_is(suppressWarnings(getNcomp(demo$block.pls, max.ncomp = 4, X=list(X=demo$X, Z=demo$Z), Y=demo$Y)), "ncomp.tune.silhouette")
    expect_is(suppressWarnings(getNcomp(demo$block.pls, max.ncomp = 4, X=list(X=demo$X, Z=demo$Z, Y = demo$Y), indY=3, scale = TRUE)), "ncomp.tune.silhouette")
})
    

test_that("getNcomp plot works",{
    res <- getNcomp(object = demo$pca, X=demo$X)
    expect_is(plot(res), "ggplot")

    res <- getNcomp(demo$pls, max.ncomp = 4, X = demo$X, Y=demo$Y)
    expect_is(plot(res), "ggplot")
    
    res <- suppressWarnings(getNcomp(demo$block.pls, max.ncomp = 4, X=list(X=demo$X, Z=demo$Z), Y=demo$Y))
    expect_is(plot(res), "ggplot")
})
