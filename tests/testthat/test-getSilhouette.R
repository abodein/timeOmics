context("getSilhouette")

demo <- suppressWarnings(get_demo_cluster())

test_that("getSilhouette failed on invalid input", {
    expect_error(getSilhouette("a"), "no applicable method for")
    expect_error(getSilhouette(1),"no applicable method for")
    expect_error(getSilhouette(c(1,2)), "no applicable method for")
    expect_error(getSilhouette(NA), "no applicable method for")
    expect_error(getSilhouette(NULL), "no applicable method for")
})

test_that("getSilhouette works", {
    expect_is(getSilhouette(demo$pca), "numeric")
    expect_is(getSilhouette(demo$spca), "numeric")
    
    expect_is(getSilhouette(demo$pls), "numeric")
    expect_is(getSilhouette(demo$spls), "numeric")
    
    expect_is(getSilhouette(demo$block.pls), "numeric")
    expect_is(getSilhouette(demo$block.spls), "numeric")
})
