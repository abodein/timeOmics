#' rebuild
#'
#' A more detailed description.
#'
#' @param X A numeric matrix.
#'
#' @param (s)pca, (s)pls, block.(s)pls results
#' that contains clustering information with molecule and cluster
#'
#'
#' @return
#' \describe{
#'   \item{One}{First item}
#'   \item{Two}{Second item}
#' }
#'
#' @examples
#'
#'
#' @export
rebuild <- function(X) UseMethod("rebuild")

rebuild.pca <- function(X){
    Y <- X$X

    # unscaling
    if(X$scale[1]){ # can be FASLE
        # warning [1] condition length > 1  => only first element will be used
        Y <- Y * X$scale
        attr(Y, "scaled:scale") <- NULL
    }

    # uncentering
    if(X$center[1]){ # can be FALSE
        # warning [1] condition length > 1  => only first element will be used
        Y.center <- matrix(rep(X$center, each = nrow(Y)), nrow = nrow(Y))
        Y <- Y + Y.center
        attr(Y, "scaled:center") <- NULL
    }
    return(Y)
}


