context("getCluster")

test_that("getCluster faild on invalid input", {
    # string
    expect_error(getCluster("abc"))
    
    # integer
    expect_error(getCluster(1))
    
    # matrix
    expect_error(getCluster(matrix(1:9, ncol = 3)))
})

demo <- suppressWarnings(get_demo_cluster())


test_that("getCluster works and return a valid output", {
    # for(i in c("pca", "spca", "pls", "spls", "block.pls", "block.spls")){
    #     expect_is(getCluster(demo[[i]]), "cluster.df")
    #     #expect_is(getCluster(demo$pca), "cluster.df")
    #}
    expect_is(getCluster(demo$pca), "cluster.df")
    expect_is(getCluster(demo$spca), "cluster.df")
    expect_is(getCluster(demo$pls), "cluster.df")
    expect_is(getCluster(demo$spls), "cluster.df")
    expect_is(getCluster(demo$block.pls), "cluster.df")
    expect_is(getCluster(demo$block.spls), "cluster.df")
    expect_is(getCluster(demo$UpDown), "cluster.df")
    
})
