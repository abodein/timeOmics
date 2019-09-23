context("test-silhouette")

test_that("params are correct", {
    data <- get_demo_silhouette()
    expect_error(silhouette(X = data$data, cluster = "z"),
                 "is(cluster, \"data.frame\") is not TRUE", fixed = TRUE)





})
