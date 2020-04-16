context("plotLong")

demo <- suppressWarnings(get_demo_cluster())
demo$block.pls2 <- suppressWarnings(mixOmics::block.pls(X=list(X=demo$X, 
                                                               Z=demo$Z), 
                                                        Y=demo$Y, ncomp = 5, 
                                                        mode = "canonical"))


pdf(NULL)
test_that("plotLong failed on invalid - object", {
    expect_error(plotLong(""), "invalid object, should be one of c(pca, spca, mixo_pls, mixo_spls, block.pls, block.spls)", fixed = TRUE)
    expect_error(plotLong(3), "invalid object, should be one of c(pca, spca, mixo_pls, mixo_spls, block.pls, block.spls)", fixed = TRUE)
    expect_error(plotLong(matrix()), "invalid object, should be one of c(pca, spca, mixo_pls, mixo_spls, block.pls, block.spls)", fixed = TRUE)    
    expect_error(plotLong(list()), "invalid object, should be one of c(pca, spca, mixo_pls, mixo_spls, block.pls, block.spls)", fixed = TRUE)
    expect_error(plotLong(data.frame()), "invalid object, should be one of c(pca, spca, mixo_pls, mixo_spls, block.pls, block.spls)", fixed = TRUE)
})

test_that("plotLong failed on invalid - time", {
    # pca / spca
    expect_error(plotLong(object = demo$pca, time = ""), "'time' should be a numeric vector")
    expect_error(plotLong(object = demo$pca, time = 1), "'time' should be a numeric vector")
    expect_error(plotLong(object = demo$pca, time = matrix()), "'time' should be a numeric vector")
    expect_error(plotLong(object = demo$pca, time = list(1:10)), "'time' should be a numeric vector")
    expect_error(plotLong(object = demo$pca, time = data.frame(1:10)), "'time' should be a numeric vector")
    expect_error(plotLong(object = demo$pca, time = NA), "'time' should be a numeric vector")
    expect_error(plotLong(object = demo$pca, time = 1:9), "'time' should be a numeric vector")
    expect_error(plotLong(object = demo$pca, time = 1:11), "'time' should be a numeric vector")

    # pls / spls
    expect_error(plotLong(object = demo$pls, time = ""), "'time' should be a numeric vector")
    expect_error(plotLong(object = demo$pls, time = 1), "'time' should be a numeric vector")
    expect_error(plotLong(object = demo$pls, time = matrix()), "'time' should be a numeric vector")
    expect_error(plotLong(object = demo$pls, time = list(1:10)), "'time' should be a numeric vector")
    expect_error(plotLong(object = demo$pls, time = data.frame(1:10)), "'time' should be a numeric vector")
    expect_error(plotLong(object = demo$pls, time = NA), "'time' should be a numeric vector")
    expect_error(plotLong(object = demo$pls, time = 1:9), "'time' should be a numeric vector")
    expect_error(plotLong(object = demo$pls, time = 1:11), "'time' should be a numeric vector")
    
    # block.pls / block.spls
    expect_error(plotLong(object = demo$block.pls, time = ""), "'time' should be a numeric vector")
    expect_error(plotLong(object = demo$block.pls, time = 1), "'time' should be a numeric vector")
    expect_error(plotLong(object = demo$block.pls, time = matrix()), "'time' should be a numeric vector")
    expect_error(plotLong(object = demo$block.pls, time = list(1:10)), "'time' should be a numeric vector")
    expect_error(plotLong(object = demo$block.pls, time = data.frame(1:10)), "'time' should be a numeric vector")
    expect_error(plotLong(object = demo$block.pls, time = NA), "'time' should be a numeric vector")
    expect_error(plotLong(object = demo$block.pls, time = 1:9), "'time' should be a numeric vector")
    expect_error(plotLong(object = demo$block.pls, time = 1:11), "'time' should be a numeric vector")
})

test_that("plotLong works", {
    expect_is(plotLong(object = demo$pca), "data.frame")
    expect_is(plotLong(object = demo$spca), "data.frame")
    expect_is(plotLong(object = demo$pls), "data.frame")
    expect_is(plotLong(object = demo$spls), "data.frame")
    expect_is(plotLong(object = demo$block.pls), "data.frame")
    expect_is(plotLong(object = demo$block.spls), "data.frame")
    expect_is(plotLong(object = demo$block.pls, time=1:10), "data.frame")
    
    # plot TRUE/FALSE
    expect_is(plotLong(object = demo$pca, plot=""), "data.frame")
    expect_is(plotLong(object = demo$pca, plot=FALSE), "data.frame")
    expect_is(plotLong(object = demo$pca, plot=1), "data.frame")
    expect_is(plotLong(object = demo$pca, plot=NULL), "data.frame")
    
    # scale TRUE/FALSE
    expect_is(plotLong(object = demo$pca, scale=""), "data.frame")
    expect_is(plotLong(object = demo$pca, scale=FALSE), "data.frame")
    expect_is(plotLong(object = demo$pca, scale=1), "data.frame")
    expect_is(plotLong(object = demo$pca, scale=NULL), "data.frame")
    
    # scale TRUE/FALSE
    expect_is(plotLong(object = demo$pca, center=""), "data.frame")
    expect_is(plotLong(object = demo$pca, center=FALSE), "data.frame")
    expect_is(plotLong(object = demo$pca, center=1), "data.frame")
    expect_is(plotLong(object = demo$pca, center=NULL), "data.frame")
    
    # legend TRUE/FALSE
    expect_is(plotLong(object = demo$pca, legend=""), "data.frame")
    expect_is(plotLong(object = demo$pca, legend=FALSE), "data.frame")
    expect_is(plotLong(object = demo$pca, legend=1), "data.frame")
    expect_is(plotLong(object = demo$pca, legend=NULL), "data.frame")
    
    # title string
    expect_is(plotLong(object = demo$pca, title=""), "data.frame")
    expect_is(plotLong(object = demo$pca, title=FALSE), "data.frame")
    expect_is(plotLong(object = demo$pca, title=1), "data.frame")
    expect_is(plotLong(object = demo$pca, title=list()), "data.frame")
    expect_is(plotLong(object = demo$pca, title=data.frame()), "data.frame")
    
    # X.label
    expect_is(plotLong(object = demo$pca, X.label=""), "data.frame")
    expect_is(plotLong(object = demo$pca, X.label=FALSE), "data.frame")
    expect_is(plotLong(object = demo$pca, X.label=1), "data.frame")
    expect_is(plotLong(object = demo$pca, X.label=list()), "data.frame")
    expect_is(plotLong(object = demo$pca, X.label=data.frame()), "data.frame")
    
    # Y.label
    expect_is(plotLong(object = demo$pca, Y.label=""), "data.frame")
    expect_is(plotLong(object = demo$pca, Y.label=FALSE), "data.frame")
    expect_is(plotLong(object = demo$pca, Y.label=1), "data.frame")
    expect_is(plotLong(object = demo$pca, Y.label=list()), "data.frame")
    expect_is(plotLong(object = demo$pca, Y.label=data.frame()), "data.frame")
    
    # legend.title
    expect_is(plotLong(object = demo$pca, legend.title=""), "data.frame")
    expect_is(plotLong(object = demo$pca, legend.title=FALSE), "data.frame")
    expect_is(plotLong(object = demo$pca, legend.title=1), "data.frame")
    expect_is(plotLong(object = demo$pca, legend.title=list()), "data.frame")
    expect_is(plotLong(object = demo$pca, legend.title=data.frame()), "data.frame")
    
    expect_is(plotLong(object = demo$pca, legend=TRUE, legend.title=""), "data.frame")
    expect_is(plotLong(object = demo$pca, legend=TRUE, legend.title=FALSE), "data.frame")
    expect_is(plotLong(object = demo$pca, legend=TRUE, legend.title=1), "data.frame")
    expect_is(plotLong(object = demo$pca, legend=TRUE, legend.title=list()), "data.frame")
    expect_is(plotLong(object = demo$pca, legend=TRUE, legend.title=data.frame()), "data.frame")
    
})
dev.off()