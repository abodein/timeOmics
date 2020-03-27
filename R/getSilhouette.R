#' Get Silhouette Coefficient from (s)PCA, (s)PLS or block.(s)PLS clusters
#'
#' \code{getSilhouette} is a generic function that compute silhouette coefficient
#' for an object of the type \code{pca}, \code{spca}, \code{pls}, \code{spls},
#' \code{block.pls}, \code{block.spls}.
#'
#' @param object a mixOmics object of the class \code{pca}, \code{spca}, \code{pls}, \code{spls}, \code{block.pls}, \code{block.spls}
#'
#' @details
#' This method extract the componant contribution depending on the object,
#' perform the clustering step, and compute the silhouette coefficient.
#'
#' @return
#' silhouette coefficient
#' 
#' @examples
#' demo <- get_demo_cluster()
#' getSilhouette(object = demo$pca)
#' getSilhouette(object = demo$spca)
#' getSilhouette(object = demo$pls)
#' getSilhouette(object = demo$spls)
#' getSilhouette(object = demo$block.pls)
#' getSilhouette(object = demo$block.spls)
#' 
#'
#' @export
getSilhouette <- function(object){
    UseMethod("getSilhouette")
}

#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @export
getSilhouette.pca <- function(object){
    cluster <- getCluster(object) %>%
        dplyr::select("molecule", "cluster")

    silhouette.res <- wrapper.silhouette(X = as.data.frame(object$X), cluster = cluster$cluster)
    return(silhouette.res$average)
}

#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @export
getSilhouette.spca <- function(object){
    cluster <- getCluster(object) %>%
        dplyr::select("molecule", "cluster")

    X <- object$X %>% as.data.frame() %>%
        dplyr::select(cluster$molecule)
    silhouette.res <- wrapper.silhouette(X = X, cluster = cluster$cluster)
    return(silhouette.res$average)
}

#' @export
#' @importFrom dplyr select
#' @importFrom magrittr %>%
getSilhouette.mixo_pls <- function(object){
    cluster <- getCluster(object) %>%
        dplyr::select("molecule", "cluster")

    X <- object$X %>% as.data.frame()
    Y <- object$Y %>% as.data.frame()
    data <- cbind(X,Y)
    silhouette.res <- wrapper.silhouette(X = data, cluster = cluster$cluster)
    return(silhouette.res$average)
}

#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @export
getSilhouette.mixo_spls <- function(object){
    cluster <- getCluster(object) %>%
        dplyr::select("molecule", "cluster")

    X <- object$X %>% as.data.frame()
    Y <- object$Y %>% as.data.frame()
    data <- cbind(X,Y) %>%
        dplyr::select(cluster$molecule)

    silhouette.res <- wrapper.silhouette(X = data, cluster = cluster$cluster)
    return(silhouette.res$average)
}

#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @export
getSilhouette.block.pls <- function(object){
    cluster <- getCluster(object) %>%
        dplyr::select("molecule", "cluster")

    X <- do.call("cbind", object$X) %>% as.data.frame

    silhouette.res <- wrapper.silhouette(X = X, cluster = cluster$cluster)
    return(silhouette.res$average)
}


#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @export
getSilhouette.block.spls <- function(object){
    cluster <- getCluster(object) %>%
        dplyr::select("molecule", "cluster")

    X <- do.call("cbind", object$X) %>% as.data.frame %>%
        dplyr::select(cluster$molecule)

    silhouette.res <- wrapper.silhouette(X = X, cluster = cluster$cluster)
    return(silhouette.res$average)
}