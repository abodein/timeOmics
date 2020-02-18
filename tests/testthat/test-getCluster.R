context("getCluster")

test_that("getCluster faild on invalid input", {
    # string
    expect_error(getCluster("abc"))
    
    # integer
    expect_error(getCluster(1))
    
    # matrix
    expect_error(getCluster(matrix(1:9, ncol = 3)))
})

demo <- get_demo_cluster()


test_that("getCluster works and return a valid output", {
    for(i in c("pca", "spca", "pls", "spls", "block.pls", "block.spls")){
        expect_is(getCluster(demo[[i]]), "data.frame")
    }
})
