context("proportionality")


test_that("proportionality faild on invalid input", {
    # string
    expect_error(proprortionality("abc"))
    
    # integer
    expect_error(proprortionality(1))
    
    # matrix
    expect_error(proprortionality(matrix(1:9, ncol = 3)))
})
    
test_that("proportionality works and return a good object", {
    demo <- suppressWarnings(suppressMessages(get_demo_cluster()))
    
    X <- list(pca=demo$pca, spca=demo$spca,
              pls=demo$pls, spls=demo$spls,
              block.pls=demo$block.pls, block.spls=demo$block.spls
    )
    for(i in X){
        expect_is(proportionality(i), "proportionality")
    }
})

test_that("proportionality plot works", {
    demo <- suppressWarnings(suppressMessages(get_demo_cluster()))
    X <- list(pca=demo$pca, spca=demo$spca,
              pls=demo$pls, spls=demo$spls,
              block.pls=demo$block.pls, block.spls=demo$block.spls)
    for(i in X){
        proportionality.res <- proportionality(i)
        plot.res <- plot(proportionality.res)
        expect_is(plot.res, "ggplot")
    }
})
