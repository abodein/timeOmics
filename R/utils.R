is_almostInteger <- function(X){
    if(!is.numeric(X) & !is.vector(X)) return(FALSE)
    if(length(X) != 1) return(FALSE)
    if(!is.finite(X)) return(FALSE)
    X.round <- round(X)
    if(X == X.round) return(TRUE)
    return(FALSE)
}

is_almostInteger_vector <- function(X){
    if(!is.vector(X) || is.list(X)){
        return(FALSE)
    }
    # if(!is.numeric(X) & !is.vector(X)) return(FALSE)
    #return(all(sapply(X, is_almostInteger)))
    return(all(vapply(X, is_almostInteger, logical(1))))
}

is_almostInteger_list <- function(X){
    if(!is.list(X)) return(FALSE)
    # if(!is.numeric(X) & !is.vector(X)) return(FALSE)
    # return(all(sapply(X, is_almostInteger)))
    return(all(vapply(X, is_almostInteger, logical(1))))
}



check_matrix <- function(X){
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


validate_matrix_X <- function(X){
    # X should be a numeric matrix
    X <- check_matrix(X)
    if(!is.numeric(X)){
        stop("X must be a numeric matrix/data.frame")
    }
    # if(any(!X)) stop("X must be a numeric matrix/data.frame")
    return(X)
}

validate_matrix_Y <- function(Y){
    # X should be a numeric matrix
    Y <- check_matrix(Y)
    if(!is.numeric(Y)){
        stop("Y must be a numeric matrix/data.frame")
    }
    # if(any(!Y)) stop("Y must be a numeric matrix/data.frame")
    return(Y)
}

validate_list_matrix_X <- function(X){
    if(!is.list(X)){
        stop("X must be a list of matrix/data.frame")
    }
    X <- lapply(X, validate_matrix_X)
    return(X)
}

validate_ncomp <- function(ncomp, X){
    # ncomp should be a positive non-null integer
    # lower than ncol(X)
    nrow_X <- ifelse(is.list(X), nrow(X[[1]]), nrow(X))
    ncomp.max <- min(unlist(lapply(X,function(x)ncol(x))), nrow_X)
    if(is.list(X)){
        ncomp.max <- min(ncomp.max, ncol(X))
    }
    if(!is_almostInteger(ncomp)){
        stop(paste0("'ncomp' should be an integer between 1 and ", ncomp.max))
    }
    if(ncomp > ncomp.max || ncomp==0){
        stop(paste0("'ncomp' should be an integer between 1 and ", ncomp.max))
    }
    ncomp <- round(ncomp)
    return(ncomp)
}

validate_test_keepX <- function(test.keepX, X){
    # test.keepX should be a vecter of positive integer of size > 1
    # every value of test.keepX should be lower than ncol(X)
    # ncomp and X have already been validate
    if(is.null(test.keepX)){
        test.keepX <- ncol(X)
    }
    if(!is_almostInteger_vector(test.keepX)){
        stop("'test.keepX' should be numeric")
    }
    if(any(test.keepX>ncol(X))){
        stop(paste0("'test.keepX' must be lower than ", ncol(X), ", ncol(X)"))
    }
    return(sort(unique(test.keepX)))
}

validate_test_keepY <- function(test.keepY, Y){
    # test.keepX should be a vecter of positive integer of size > 1
    # every value of test.keepX should be lower than ncol(X)
    # ncomp and X have already been validate
    if(is.null(Y)){  # case of keepY null in block spls
        return(NULL)
    } else {
        if(is.null(test.keepY)){
            test.keepY <- ncol(Y)
        }
        if(!is_almostInteger_vector(test.keepY)){
            stop("'test.keepY' should be numeric")
        }
        if(any(test.keepY>ncol(Y))){
            stop(paste0("'test.keepY' must be lower than ", ncol(Y), ", ncol(Y)"))
        }
    }
    return(sort(unique(test.keepY)))
}

validate_test_list_keepX <- function(test.keepX, ncomp, X){
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
    if(is_almostInteger_list(test.keepX)){
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
    test.keepX <- lapply(test.keepX, function(x){x})
    return(test.keepX)
}

validate_indY <- function(indY, X){
    # X already checked
    if(is.null(indY)){
        stop(paste0("'indY' must be a numeric value lower or equal to ", length(X), ", the number of blocks in X."))
    }
    if(!is_almostInteger(indY) | !(indY %in% c(1:length(X))) ){
        stop(paste0("'indY' must be a numeric value lower or equal to ", length(X), ", the number of blocks in X."))
    }
    return(indY)
}

sd_new <- function(x, ...){
    if(length(x) == 1){
        return(0)
    }else{
        return(sd(x, ...))
    }
}

return_true_false <- function(x, default){
    if(is.logical(x)){
        if(is.finite(x)){
            return(x)
        } else { #NA
            return(default)
        }
    } else {
        return(default)
    }
}

check_legend.block.name <- function(legend.block.name, cluster){
    stopifnot(is(legend.block.name, "character"))
    stopifnot(is(legend.block.name, "vector"))
    
    # cluster <- getCluster(...)
    size_block <- unique(cluster$block)
    stopifnot(length(legend.block.name) == length(size_block))
}
