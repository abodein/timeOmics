# test <- mixOmics::pca(demo$X, ncomp = 2)
#
# C <- test$call
# C$ncomp <- 3
# H <- eval(C)
#
# demo <- suppressMessages(get_demo_cluster())
# spca.res <- demo$spca
# spca.cluster <- getCluster(spca.res)
# plotLong(spca.res, scale = T)
