context("getNcomp")

demo <- get_demo_cluster()
demo$pca <- mixOmics::pca(X = demo$X, ncomp = 5)
demo$pls <- mixOmics::pls(X = demo$X, Y=demo$Y, ncomp =5, mode = "canonical")
# demo$block.pls = mixOmics::block.pls(X=list(X=demo$X, Y=demo$Y, Z=demo$Z), indY=1, ncomp = 5, mode = "canonical")


test_that("getNcomp failed on invalid input", {
    # test for "object"
    # string, integer, matrix
    lapply(list("",1,demo$X), function(i) expect_error(getNcomp(i)))
    
    # test for "max.ncomp"
    # wrong input class
    lapply(list("", NA), function(i)
        expect_error(getNcomp(object = demo$pca, max.ncomp = i, ncomp=5), "invalid value for 'max.ncomp'."))
    expect_error(getNcomp(object = demo$pca, max.ncomp = demo$X, ncomp=5), "length(max.ncomp) == 1 is not TRUE", fixed = TRUE)
    # wrong value
    expect_error(getNcomp(object = demo$pca, max.ncomp = 11, ncomp=5), "use smaller 'max.ncomp'")
    
    # Missing parameters
    expect_error(getNcomp(object = demo$pca), "Missing parameters, please provide the same parameters as the ones contained in the mixOmics object.", fixed = TRUE)
})
# 
# test_that("getNcomp plot works",{
#     res <- getNcomp(object = demo$pca, X=demo$X, ncomp = 5)
#     expect_is(plot(res), "numeric")
# 
# })

# test_that("getNcomp works and return valid output class", {
#     # pca 
    # expect_is(getNcomp(object = demo$pca, X=demo$X, ncomp = 5), "ncomp.tune.silhouette")
    #res <- getNcomp(object = demo$pca, X=demo$X, ncomp = 5)
    #expect_is(res, "ncomp.tune.silhouette")
    # expect_true(is(res))
    
#    # pls
#    expect_is(getNcomp(object = demo$pls, X=demo$X, Y=demo$Y, ncomp = 5, mode = "canonical"), "ncomp.tune.silhouette")
#    # block.pls
#    expect_is(getNcomp(object = demo$block.pls, X=list(X=demo$X, Y=demo$Y, Z=demo$Z), indY=1, ncomp = 5, mode = "canonical"), "ncomp.tune.silhouette")
# })

# debug(getNcomp)
# # 
# expect_is(getNcomp(demo$pca), "ncomp.tune.silhouette")
# expect_is(getNcomp(demo$pls, X= demo$X, Y = demo$Y, mode ="canonical", ncomp = 5), "ncomp.tune.silhouette")
# expect_is(getNcomp(demo$block.pls), "ncomp.tune.silhouette")
