#' Remove features with low variation
#'
#'
#' \code{remove.low.cv} that removes variables with low variation.
#' From a matrix/data.frame (samples in rows, features in columns), it computes the coefficient of variation for every features (columns) 
#' and return a filtered data.frame with features for which the coefficient of variation is above a given threshold.
#'
#' @param X a matrix/data.frame
#' @param cutoff a numeric value
#'
#' @return
#' a data.frame/matrix
#' 
#' @examples 
#' mat <- matrix(sample(1:3, size = 200, replace = TRUE), ncol=20)
#' remove.low.cv(mat, 0.4)
#' 
#' @export
remove.low.cv <- function(X, cutoff = 0.5){
    stopifnot(is(X, "data.frame") | is(X, "matrix"))
    stopifnot(is.vector(cutoff) & is.numeric(cutoff) & (length(cutoff) == 1))
    
    # var.coef
    cv <- unlist(lapply(as.data.frame(X), 
                        function(x) abs(sd(x)/mean(x))))
    return(X[,cv > cutoff])
}
