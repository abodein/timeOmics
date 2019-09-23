#' Get data for silhouette demo
#'
#' @return A matrix of expression profile, sample in raws, time in columns.
#'
#' @examples
#' data <- get_demo_silhouette()
#'
#' @export
get_demo_silhouette <- function() {
    readRDS(system.file("extdata/data_silhouette.rds", package="timeOmics",
                        mustWork = TRUE))
}


#' silhouette
#'
#' A more detailed description.
#'
#' @param X A numeric matrix of size NxP with feature in colnames.
#'
#' @param cluster A numeric vector of sive ncol(X) containing the cluster
#' information by
#'
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @export
silhouette <- function(dmatrix,  # distance matrix
                       cluster)  # cluster vector of size ncol(dmatrix)
    {
    #-- checking input parameters ---------------------------------------------#
    #--------------------------------------------------------------------------#

    #-- check dmatrix
    stopifnot(is.matrix(dmatrix) || is.data.frame(dmatrix))
    dmatrix <- as.matrix(dmatrix)
    stopifnot(nrow(dmatrix) == ncol(dmatrix))

    #-- check cluster
    stopifnot(is.vector(cluster) || is.factor(cluster))
    stopifnot(length(cluster) == ncol(dmatrix))

    cluster <- factor(cluster)
    cluster.levels <- levels(cluster)


    #- compute silhouette -----------------------------------------------------#
    #--------------------------------------------------------------------------#
    average.dist <- matrix(ncol = length(cluster.levels), nrow = nrow(dmatrix))
    result <- vector(length = ncol(dmatrix))
    for(i in 1:nrow(dmatrix)){
        for(j in 1:length(cluster.levels)){
            index.tmp <- cluster == cluster.levels[j]
            if(cluster.levels[j] == cluster[i]) {
                # we do not include d(i,i) in the sum
                # can introduce NaN if size of cluster is 1
                # but this is handle after because score is 0 when size of cluster is 1
                index.tmp[i] <- FALSE
            }
            average.dist[i,j] <- mean(dmatrix[i, index.tmp])
        }
        A <- average.dist[i, cluster[i]] # a : inside
        B <- min(average.dist[i,-c(cluster[i])])  # b
        result[i] <- silhoutte.formula(A = A, B = B)
    }

    #-- return
    to_return <- list()

    #-- silhouette coefficient by feature
    to_return[["feature"]] <- cbind(colnames(dmatrix), cluster,  as.data.frame(result))
    colnames(to_return[["feature"]]) <- c("feature", "cluster", "silhouette.coef")

    #-- average silhouette coefficient
    to_return[["average"]] <- mean(to_return[["feature"]][["silhouette.coef"]])

    #-- average silhouette coefficient by cluster
    to_return[["average.cluster"]] <- dplyr::group_by(to_return[["feature"]], cluster) %>%
        dplyr::summarise(silhouette.coef = mean(silhouette.coef)) %>%
        as.data.frame

    return(to_return)
}

#' dmatrix.spearman.dissimilarity
#'
#' Compute the spearman dissimilarity distance.
#'
#' @param X A numeric matrix with feature in colnames
#'
#' @return
#' Return a dissimilarity matrix of size PxP.
#'
#' @export
dmatrix.spearman.dissimilarity <- function(X){
    # between 0 and 2
    dmatrix <- cor(x = X, use = 'pairwise.complete.obs',
                   method = 'spearman')
    dmatrix <- 1 - dmatrix
    return(dmatrix)
}

dmatrix.proportionality.distance <- function(X){
    # clr first
    dmatrix <- matrix(ncol = ncol(X), nrow = ncol(X))
    rownames(dmatrix) <- colnames(dmatrix) <- colnames(X)
    for(i in 1:ncol(X)){
        for(j in 1:ncol(X)){
            dmatrix[i,j] <- var(X[,i] - X[,j])/var(X[,i] + X[,j])
        }
    }
    return(dmatrix)
}

silhoutte.formula <- function(A,B){
    # A average dist inside cluster; # B: min average dist outside
    stopifnot(is.vector(A) && is.vector(B))
    stopifnot(is.numeric(A) && is.numeric(B))
    if(!(is.finite(A))){
        return(0)
    }else{
        return((B - A)/(max(A,B)))
    }
}

#' @export
wrapper.silhouette <- function(X, cluster)
{
    #-- checking input parameters ---------------------------------------------#
    #--------------------------------------------------------------------------#

    #-- check X
    stopifnot(is.matrix(X) || is.data.frame(X))
    X <- as.matrix(X)

    #-- check cluster
    stopifnot(is.vector(cluster) || is.factor(cluster))
    stopifnot(length(cluster) == ncol(X))

    cluster <- factor(cluster)
    cluster.levels <- levels(cluster)

    #-- compute distance matrix -----------------------------------------------#
    #--------------------------------------------------------------------------#
    dmatrix <- dmatrix.spearman.dissimilarity(X)

    #- compute silhouette -----------------------------------------------------#
    #--------------------------------------------------------------------------#
    silhouette(dmatrix = dmatrix, cluster = cluster)
}
