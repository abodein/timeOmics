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
    check.matrix(X)
    return(X)
}

validate.list.matrix.X <- function(X){
    X <- lapply(X, validate.matrix.X)
    return(X)
}

validate.ncomp <- function(ncomp, X){
    # ncomp should be a positive non-null integer
    # lower than ncol(X) - 1
    ncomp <- round(ncomp)
    return(ncomp)
}

validate.test.keepX <- function(test.keepX, ncomp, X){
    # test.keepX should be a vecter of positive integer of size > 1
    # every value of test.keepX should be lower than ncol(X)
    # ncomp and X have already been validate
    if(is.null(X)){  # case of keepY null in block spls
        return(NULL)
    }
    return(test.keepX)
}

validate.test.list.keepX <- function(test.keepX, ncomp, X){
    # for block spls
    # same length of X (list)
    return(test.keepX)
}

validate.indY <- function(indY, X){
    if(!is.almostInteger(indY)){
        return(FALSE)
    }
    if(!(indY %in% c(1:length(X)))){
        # data is valid and is a list of df 
        return(FALSE)
    }
    return(indY)
}


