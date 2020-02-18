#' Get optimal number of components
#'
#' Compute the average silhouette coefficient for a given set of components on a mixOmics result.
#' Foreach given ncomp, the mixOmics method is performed with the sames arguments and the given `ncomp`.
#' Longitudinal clustering is performed and average silhouette coefficient is computed.
#'
#' @param object A mixOmics object of the class `pca`, `spca`, `mixo_pls`, `mixo_spls`, `block.pls`, `block.spls`
#'
#' @param min.ncomp integer, minimum number of component to include.
#' Must be greater than 0. If no argument is given, `min.ncomp = 1`.
#'
#' @param max.ncomp integer, minimum number of component to include.
#' If no argument is given, `max.ncomp=object$ncomp`
#'
#' @return
#' \code{getNcomp} returns a list with class "ncomp.tune.silhouette" containing the following components:
#'
#' \item{ncomp}{a vector containing the tested ncomp}
#' \item{silhouette}{a vector containing the average silhouette coefficient by ncomp}
#' \item{dmatrix}{the distance matrix used to compute silhouette coefficient}
#'
#' @seealso \code{\link{getCluster}}, \code{\link{silhouette}}
#'
#' @examples
#' # random input data
#' X <- matrix(rnorm(1:100), ncol = 10, dimnames = list(1:10, paste0("X_", 1:10)))
#'
#' # Principal Component Analysis
#' pca.res <- pca(X, ncomp = 5)
#'
#' # getNcomp
#' res.ncomp <- getNcomp(pca.res, max.ncomp = 4)
#' plot(res.ncomp)
#'
#' @export
#'
#' @import mixOmics
getNcomp <- function(object, max.ncomp = NULL, X, Y = NULL, ...){
    #-- checking input parameters ---------------------------------------------#
    #--------------------------------------------------------------------------#

    #-- object
    stopifnot( is(object, c("pca", "mixo_pls", "block.pls")))

    #-- max.ncomp
    if(!is.null(max.ncomp)){
        stopifnot(length(max.ncomp) == 1)
        if ( !is.numeric(max.ncomp) || max.ncomp < 1 || !is.finite(max.ncomp))
            stop("invalid value for 'max.ncomp'.")

        if (max.ncomp > min(ncol(object$X), nrow(object$X)))
            stop("use smaller 'max.ncomp'")
    } else {
        max.ncomp <- unique(object$ncomp)
    }

    #-- check for correct parameters in object$call
    Args = as.list(match.call())
    mixo.call <- object$call
    #mixo.call[[1]] <- as.name(paste0("mixOmics::", mixo.call[[1]]))
    if(!(all(names(mixo.call)[-1] %in% names(Args)[-1]))){
        stop("Missing parameters, please provide the same parameters as the ones contained in the mixOmics object.")
    }
    
    #-- Iterating ncomp
    silhouette.res <- vector(length = max.ncomp)

    #-- compute dmatrix using spearman dissimilarity
    XX <- object$X
    if(is.null(dim(XX))){
        XX <- do.call("cbind", XX)
    }
    # if Y
    if(!is.null(object$Y)){
        XX <- cbind(XX,object$Y)
    }
    
    dmatrix <- dmatrix.spearman.dissimilarity(XX)

    #-- iterative
    i <- 1
    for(comp in 1:max.ncomp){
        mixo.call$ncomp <- comp
        mixo.res <- eval(mixo.call)
        cluster <- getCluster(mixo.res)
        # order of feature is the same as colnames(X)
        cluster %>% mutate(molecule = factor(molecule, levels = colnames(dmatrix))) 
        stopifnot(cluster$molecule == colnames(dmatrix))
        silhouette.res[i] <- silhouette(dmatrix, cluster$cluster)$average
        i <- i + 1
    }

    to_return <- list()
    to_return[["ncomp"]] <- c(0,1:max.ncomp)
    to_return[["silhouette"]] <- c(0,silhouette.res)
    to_return[["dmatrix"]] <- dmatrix

    class(to_return) <- "ncomp.tune.silhouette"
    return(invisible(to_return))
}

#' @export
plot.ncomp.tune.silhouette <- function(object){
    stopifnot(is(object, "ncomp.tune.silhouette"))
    plot(x = object$ncomp, y = object$silhouette, type = "b", xaxt="n",
         xlab = "Number of Principal Components", ylab = "Average Silhouette Coefficient")
    axis(side = 1, at = object$ncomp)
}
