#' getSilhouette
#'
#' \code{getSilhouette} is a generic function that compute silhouette coefficient
#' for an object of the type \code{pca}, \code{spca}, \code{pls}, \code{spls},
#' \code{block.pls}, \code{block.spls}.
#'
#' @param X A numeric matrix.
#'
#' @param cluster A data.frame
#' that contains clustering information with molecule and cluster
#'
#' @details
#' This method extract the componant contribution depending on the object,
#' perform the clustering step, and compute the silhouette coefficient.
#'
#' @return
#' 
#' @examples
#' demo <- suppressMessages(get_demo_cluster())
#' # pca
#' getSilhouette(object = demo$pca)
#' # spca
#' getSilhouette(object = demo$spca)
#' 
#'
#' @export
getSilhouette <- function(object){
    UseMethod("getSilhouette")
}

#' @importFrom magrittr %>%
#' @importFrom dplyr select
getSilhouette.pca <- function(object){
    cluster <- getCluster(object) %>%
        dplyr::select("molecule", "cluster")

    silhouette.res <- wrapper.silhouette(X = as.data.frame(object$X), cluster = cluster$cluster)
    return(silhouette.res$average)
}

#' @importFrom magrittr %>%
#' @importFrom dplyr select
getSilhouette.spca <- function(object){
    cluster <- getCluster(object) %>%
        dplyr::select("molecule", "cluster")

    X <- object$X %>% as.data.frame() %>%
        dplyr::select(cluster$molecule)
    silhouette.res <- wrapper.silhouette(X = X, cluster = cluster)
    ComputeSilhouetteAverage(silhouette.res)
}

#' @importFrom magrittr %>%
getSilhouette.mixo_pls <- function(object){
    cluster <- getCluster(object) %>%
        dplyr::select("molecule", "cluster")

    X <- object$X %>% as.data.frame()
    Y <- object$Y %>% as.data.frame()
    data <- cbind(X,Y)
    silhouette.res <- silhouette(X = data, cluster = cluster)
    ComputeSilhouetteAverage(silhouette.res)
}

#' @importFrom magrittr %>%
#' @importFrom dplyr select
getSilhouette.mixo_spls <- function(object){
    cluster <- getCluster(object) %>%
        dplyr::select("molecule", "cluster")

    X <- object$X %>% as.data.frame()
    Y <- object$Y %>% as.data.frame()
    data <- cbind(X,Y) %>%
        dplyr::select(cluster$molecule)

    silhouette.res <- silhouette(X = data, cluster = cluster)
    ComputeSilhouetteAverage(silhouette.res)
}

#' @importFrom magrittr %>%
getSilhouette.block.pls <- function(object){
    cluster <- getCluster(object) %>%
        dplyr::select("molecule", "cluster")

    X <- do.call("cbind", object$X) %>% as.data.frame

    silhouette.res <- silhouette(X = X, cluster = cluster)
    ComputeSilhouetteAverage(silhouette.res)
}


#' @importFrom magrittr %>%
#' @importFrom dplyr select
getSilhouette.block.pls <- function(object){
    cluster <- getCluster(object) %>%
        dplyr::select("molecule", "cluster")

    X <- do.call("cbind", object$X) %>% as.data.frame %>%
        dplyr::select(cluster$molecule)

    silhouette.res <- silhouette(X = X, cluster = cluster)
    ComputeSilhouetteAverage(silhouette.res)
}
