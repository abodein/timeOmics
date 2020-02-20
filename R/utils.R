#' check numeric matrix
#'
is.numeric.matrix <- function(X){
    if(any(is.infinite(X))) return(FALSE)
    if(is.numeric(X) & is.matrix(X)) return(TRUE)
    return(FALSE)
}

is.almostInteger <- function(X){
    if(!is.numeric(X) & !is.vector(X)) return(FALSE)
    if(length(X) != 1) return(FALSE)
    if(!is.finite(X)) return(FALSE)
    X.round <- round(X)
    if(X == X.round) return(TRUE)
    return(FALSE)
}

is.almostInteger.vector <- function(X){
    if(is.list(X)){
        return(FALSE)
    }
    # if(!is.numeric(X) & !is.vector(X)) return(FALSE)
    return(all(sapply(X, is.almostInteger)))
}

is.almostInteger.list <- function(X){
    if(!is.list(X)) return(FALSE)
    # if(!is.numeric(X) & !is.vector(X)) return(FALSE)
    return(all(sapply(X, is.almostInteger)))
}



check.matrix <- function(X){
    # add rownames and colnames if absent, cast into matrix
    if(!(is.matrix(X) || is.data.frame(X))) return(FALSE)

    if(is.data.frame(X)){
        X <- as.matrix(X)
    }
    if(is.null(rownames(X))){
        rownames(X) <- 1:nrow(X)
    }
    if(is.null(colnames(X))){
        colnames(X) <- paste0("V", 1:ncol(X))
    }
    return(X)
}


validate.matrix.X <- function(X){
    # X should be a numeric matrix
    X <- check.matrix(X)
    suppressWarnings(if(!X) stop("X must be a numeric matrix/data.frame"))
    return(X)
}

validate.matrix.Y <- function(X){
    # X should be a numeric matrix
    X <- check.matrix(X)
    suppressWarnings(if(!X) stop("Y must be a numeric matrix/data.frame"))
    return(X)
}

validate.list.matrix.X <- function(X){
    if(!is.list(X)){
        stop("X must be a list of matrix/data.frame")
    }
    X <- lapply(X, validate.matrix.X)
    return(X)
}

validate.ncomp <- function(ncomp, X){
    # ncomp should be a positive non-null integer
    # lower than ncol(X) - 1
    ncomp.max <- min(unlist(lapply(X,function(x)ncol(x)-1)))
    if(!is.almostInteger(ncomp)){
        stop(paste0("'ncomp' should be an integer between 1 and ", ncomp.max, ", min(unlist(lapply(X,function(x)ncol(x)-1)))"))
    }
    if(ncomp > ncomp.max){
        stop(paste0("'ncomp' should be an integer between 1 and ", ncomp.max, ", min(unlist(lapply(X,function(x)ncol(x)-1)))"))
    }
    ncomp <- round(ncomp)
    return(ncomp)
}

validate.test.keepX <- function(test.keepX, ncomp, X){
    # test.keepX should be a vecter of positive integer of size > 1
    # every value of test.keepX should be lower than ncol(X)
    # ncomp and X have already been validate
    if(is.null(test.keepX)){
        stop("'test.keepX' should be numeric")
    }
    if(!is.almostInteger.vector(X))
    return(test.keepX)
}

validate.test.keepY <- function(test.keepY, Y){
    # test.keepX should be a vecter of positive integer of size > 1
    # every value of test.keepX should be lower than ncol(X)
    # ncomp and X have already been validate
    if(is.null(Y)){  # case of keepY null in block spls
        return(NULL)
    } else {
        if(is.null(test.keepY)){
            test.keepY <- ncol(Y)
        }
        if(!is.almostInteger.vector(test.keepY))
            stop("'test.keepY' should be numeric")
    }
    return(test.keepY)
}

validate.test.list.keepX <- function(test.keepX, ncomp, X){
    # for block spls
    # same length of X (list)
    # if (is.null(test.keepX)) {
    #     test.keepX = lapply(seq_along(X), function(x) {
    #         c(5, 10, 15)[which(c(5, 10, 15) < ncol(X[[x]]))]
    #     })
    #     names(test.keepX) = names(X)
    # }
    if(is.null(test.keepX)){
        stop(paste0("'test.list.keepX' must be a list of numeric of size ", length(X), "."))
    }
    if(is.almostInteger.list(test.keepX)){
        stop(paste0("'test.list.keepX' must be a list of numeric of size ", length(X), "."))
    }
    if(!(all(names(test.keepX) %in% names(X)) && all(names(X) %in% names(test.keepX)))){
        stop("'list.test.keepX' should have the same names as X")
    }
    lapply(1:length(X), function(i){
        if(any(ncol(X[[i]]) < test.keepX[[i]])){
            stop(paste0("'test.list.keepX[[",i,"]] sould be lower than ",ncol(X[[i]]),", ncol(X[[",i,"]])."))
        }
    })
    return(test.keepX)
}

validate.indY <- function(indY, X){
    # X already checked
    if(is.null(indY)){
        stop(paste0("'indY' must be a numeric value lower or equal to ", length(X), ", the number of blocks in X."))
    }
    if(!is.almostInteger(indY) | !(indY %in% c(1:length(X))) ){
        stop(paste0("'indY' must be a numeric value lower or equal to ", length(X), ", the number of blocks in X."))
    }
    return(indY)
}
