#' Get variable cluster from (s)PCA, (s)PLS or block.(s)PLS
#'
#' This function returns the cluster associated to each feature from a mixOmics object.
#'
#' @param X an object of the class: \code{pca}, \code{spca}, \code{pls}, \code{spls}, \code{block.pls} or \code{block.spls}
#'
#' @return
#' A data.frame containing the name of feature, its assigned cluster and other information such as selected component, contribution, sign, ...
#'
#' @details
#' For each feature, the cluster is assigned according to the maximum contribution on a component and the sign of that contribution.
#' 
#' @seealso 
#' \code{\link[mixOmics]{selectVar}}
#' 
#' @examples
#' demo <- suppressWarnings(get_demo_cluster())
#' pca.cluster <- getCluster(demo$pca)
#' spca.cluster <- getCluster(demo$spca)
#' pls.cluster <- getCluster(demo$pls)
#' spls.cluster <- getCluster(demo$spls)
#' block.pls.cluster <- getCluster(demo$block.pls)
#' block.spls.cluster <- getCluster(demo$block.spls)
#'
#' @export
getCluster <- function(X) UseMethod("getCluster")

#' get_demo_cluster
#' 
#' Generates random data to be used in examples.
#' 
#' @return  a list containg:
#' \item{X}{data.frame}
#' \item{Y}{data.frame}
#' \item{Z}{data.frame}
#' \item{pca}{a mixOmics pca result}
#' \item{spca}{a mixOmics spca result}
#' \item{pls}{a mixOmics pls result}
#' \item{spls}{a mixOmics spls result}
#' \item{block.pls}{a mixOmics block.pls result}
#' \item{block.spls}{a mixOmics block.spls result}
#' 
#' @examples
#' # Random data could lead to "The SGCCA algorithm did not converge" warning which is not important for a demo
#' demo <- suppressWarnings(get_demo_cluster())
#' @export
get_demo_cluster<- function(){
    X <- matrix(sample(1:1000), nrow = 10,
                dimnames = list(1:10, paste0("X_",1:100)))

    Y <- matrix(sample(1:100), nrow = 10, 
                dimnames = list(1:10, paste0("Y_",1:10)))

    Z <- matrix(sample(1:500), nrow = 10, 
                dimnames = list(1:10, Y = paste0("Z_",1:50)))
    
    list.res = list()
    list.res$X <- X
    list.res$Y <- Y
    list.res$Z <- Z
    list.res$pca <- mixOmics::pca(X = X, ncomp = 5)
    list.res$spca <- mixOmics::spca(X = X, ncomp = 5, keepX = c(5, 15, 4,5,6))

    list.res$pls <- mixOmics::pls(X = X, Y = Y, ncomp = 5, mode = "canonical")
    list.res$spls <- mixOmics::spls(X = X, Y = Y, ncomp = 5, mode = "canonical",
                                keepX = c(5,6,4,5,6), keepY = c(5,1,4,5,6))

    list.res$block.pls <- mixOmics::block.pls(X = list("X" = X, "Y" = Y, "Z" = Z), indY = 1,
                                             ncomp = 5, mode = "canonical")

    list.res$block.spls <- mixOmics::block.spls(X = list("X" = X, "Y" = Y, "Z" = Z), indY = 1, ncomp = 3,
                                             mode = "canonical", keepX = list("X" = c(5,6,4), "Y" = c(5,5,5), "Z" = c(4,2,4)))
    return(invisible(list.res))
}

#' @importFrom dplyr mutate
#' @importFrom tibble rownames_to_column
#' @importFrom stringr str_remove
#' @importFrom magrittr %>%
#' @export
getCluster.pca <- function(X){
    #print("getCluster.pca")
    # colnames = PC1, PC2...
    loadings.max <- getMaxContrib(X$loadings$X)

    loadings.max <- loadings.max %>% 
        rownames_to_column("molecule") %>%
        mutate(cluster = stringr::str_remove(comp, "^PC") %>% 
                   as.numeric()) %>%
        mutate(block = "X") %>%
        .mutate_cluster()
    Valid.getCluster(loadings.max)
    return(loadings.max)
}

#' @export
#' @importFrom dplyr mutate
#' @importFrom tibble rownames_to_column
#' @importFrom stringr str_remove
#' @importFrom magrittr %>%
getCluster.spca <- function(X){
    # print(class(X))
    selected.features.loadings <- X$loadings$X[rowSums(X$loadings$X) != 0,,drop=FALSE]
    loadings.max <- getMaxContrib(selected.features.loadings)

    loadings.max <- loadings.max %>% 
        rownames_to_column("molecule") %>%
        mutate(cluster = stringr::str_remove(comp, "^PC") %>% 
                   as.numeric()) %>%
        mutate(block = "X") %>%
        .mutate_cluster()
    Valid.getCluster(loadings.max)
    return(loadings.max)
}

#' @export
#' @importFrom dplyr mutate
#' @importFrom tibble rownames_to_column
#' @importFrom stringr str_remove
#' @importFrom magrittr %>%
getCluster.mixo_pls <- function(X){
    # print(class(X))
    # block X
    loadings.max.X <- getMaxContrib(X$loadings$X)

    loadings.max.X <- loadings.max.X %>% 
        rownames_to_column("molecule") %>%
        mutate(cluster = stringr::str_remove(comp, "^comp") %>% 
                   as.numeric()) %>%
        mutate(block = "X")


    # block Y
    loadings.max.Y <- getMaxContrib(X$loadings$Y)

    loadings.max.Y <- loadings.max.Y %>% 
        rownames_to_column("molecule") %>%
        mutate(cluster = stringr::str_remove(comp, "^comp") %>% 
                   as.numeric()) %>%
        mutate(block = "Y")

    loadings.max <- rbind(loadings.max.X, loadings.max.Y) %>%
        .mutate_cluster()

    Valid.getCluster(loadings.max)
    return(loadings.max)
}

#' @export
#' @importFrom dplyr mutate
#' @importFrom tibble rownames_to_column
#' @importFrom stringr str_remove
#' @importFrom magrittr %>%
getCluster.mixo_spls <- function(X){
    # note : can not concatenate X and Y
    # because they can have the same features names contrary to block.(s)pls

    # print(class(X))
    # block X
    X.selected.features.loadings <- X$loadings$X[rowSums(X$loadings$X) != 0,,drop=FALSE]
    loadings.max.X <- getMaxContrib(X.selected.features.loadings)

    loadings.max.X <- loadings.max.X %>%
        rownames_to_column("molecule") %>%
        mutate(cluster = stringr::str_remove(comp, "^comp") %>%
                   as.numeric()) %>%
        mutate(block = "X")

    # block Y
    Y.selected.features.loadings <- X$loadings$Y[rowSums(X$loadings$Y) != 0,]
    loadings.max.Y <- getMaxContrib(Y.selected.features.loadings)

    loadings.max.Y <- loadings.max.Y %>%
        rownames_to_column("molecule") %>%
        mutate(cluster = stringr::str_remove(comp, "^comp") %>%
                   as.numeric()) %>%
        mutate(block = "Y")

    loadings.max <- rbind(loadings.max.X, loadings.max.Y)  %>%
        .mutate_cluster()

    Valid.getCluster(loadings.max)
    return(loadings.max)
}

#' @export
#' @importFrom purrr imap set_names
#' @importFrom dplyr mutate left_join
#' @importFrom tibble rownames_to_column
#' @importFrom stringr str_remove
#' @importFrom magrittr %>%
getCluster.block.pls <- function(X){
    # print(class(X))
    # get block info
    block.info <- purrr::imap(X$loadings, function(x,y) rownames(x) %>%
                           as.data.frame %>%
                           set_names("molecule") %>%
                           mutate("block" = y))
    block.info <- do.call("rbind", block.info) %>% as.data.frame() %>%
        mutate(block = factor(block, levels = names(X$loadings))) %>%
        mutate(molecule = as.character(molecule))

    loadings <- do.call("rbind", X$loadings)
    loadings.max <- getMaxContrib(loadings)

    loadings.max <- loadings.max %>%
        rownames_to_column("molecule") %>%
        mutate(cluster = stringr::str_remove(comp, "^comp") %>%
                   as.numeric()) %>%
        left_join(block.info, by = c("molecule", "molecule")) %>%
        .mutate_cluster()

    Valid.getCluster(loadings.max)
    return(loadings.max)
}

#' @export
#' @importFrom purrr imap set_names
#' @importFrom dplyr mutate left_join
#' @importFrom tibble rownames_to_column
#' @importFrom stringr str_remove
#' @importFrom magrittr %>%
getCluster.block.spls <- function(X){

    # print(class(X))
    # get block info
    block.info <- purrr::imap(X$loadings, function(x,y) rownames(x) %>%
                                  as.data.frame %>%
                                  set_names("molecule") %>%
                                  mutate("block" = y))
    block.info <- do.call("rbind", block.info) %>%
        as.data.frame() %>%
        mutate(block = factor(block, levels = names(X$loadings))) %>%
        mutate(molecule = as.character(molecule))

    # sparse
    loadings <- do.call("rbind", X$loadings)
    X.selected.features.loadings <- loadings[rowSums(loadings) != 0,, drop=FALSE]
    loadings.max <- getMaxContrib(X.selected.features.loadings)

    loadings.max <- loadings.max %>%
        rownames_to_column("molecule") %>%
        mutate(cluster = stringr::str_remove(comp, "^comp") %>%
                   as.numeric()) %>%
        left_join(block.info, by = c("molecule", "molecule")) %>%
        .mutate_cluster()

    Valid.getCluster(loadings.max)
    return(loadings.max)
}

#' Get Max Contrib from loading matrix
#' 
#' @param X loading matrix from mixOmics
#' @return a matrix
#' 
#' @keywords internal
#' @noRd 
#' @importFrom purrr set_names
#' @importFrom magrittr %>%
getMaxContrib <- function(X){
    # loadings matrix, features in rows, comp in columns
    contrib.max <- apply(X = X, FUN = function(x) { x[which.max( abs(x) )][1]}, MARGIN = 1) %>%
        as.data.frame() %>%
        purrr::set_names("contrib.max")

    cluster.info <- apply(X = X, FUN = function(x) { colnames(X)[which.max( abs(x) )[1]]}, MARGIN = 1) %>%
        as.data.frame() %>%
        purrr::set_names("comp")

    stopifnot(rownames(contrib.max) == rownames(cluster.info))
    return(cbind(cluster.info, contrib.max))
}

# absmax <- function(x) { x[which.max( abs(x) )][1]}
# absmax.index <- function(x) { which.max( abs(x) )[1]}


#' @importFrom dplyr mutate case_when pull
#' @importFrom magrittr %>%
.mutate_cluster <- function(loadings.max){
    X <- loadings.max %>%
        mutate(cluster = cluster * sign(contrib.max)) %>%
    mutate(contribution = case_when(sign(contrib.max) == 1 ~ "positive",
                                    sign(contrib.max) == -1 ~ "negative",
                                    sign(contrib.max) == 0 ~ "NULL"))

    cluster.order <-  X %>%
        pull(cluster) %>%
        abs %>%
        unique %>% 
        sort %>%
        rep(each=2) %>% `*`(c(1,-1))

    X <- X %>%
        mutate(cluster = factor(cluster, levels = cluster.order))
    return(X)
}

Valid.getCluster <- function(X){
    col_names <- c("molecule","comp","contrib.max","cluster","block","contribution")
    stopifnot(all(col_names %in% colnames(X)))
    # other check ?
    # all comp present? sometimes not true
    # idem for number of cluster
    # also a molecule can be found in different cluster
}
