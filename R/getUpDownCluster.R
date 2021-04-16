#' Up-Down clustering
#' 
#' Performs a clustering based on the signs of variation between 2 timepoints.
#' 
#' @param X a dataframe or list of dataframe with the same number of rows.
#' 
#' @examples
#' demo <- suppressWarnings(get_demo_cluster())
#' X <- list(X = demo$X, Y = demo$Y, Z = demo$Z)
#' res <- getUpDownCluster(X)
#' class(res)
#' getCluster(res)
#' 
#' X <- demo$X
#' res <- getUpDownCluster(X)

#' @importFrom purrr imap_dfr
#' @export
getUpDownCluster <- function(X){

    #stopifnot(class(X) %in% c("matrix", "data.frame", "list"))
    stopifnot(is(X, "matrix") || is(X, "data.frame") || is(X, "list"))
    
    
    
    if(is.matrix(X) || is.data.frame(X)){
        
        # check X
        X <- validate_matrix_X(X)
        X <- as.data.frame(X)
        
        res <- .getUpDown(X) %>% mutate(block = "X")
    }
    else if(is.list(X) & length(X)>1){
        
        # check X list
        X <- validate_list_matrix_X(X)
        X <- lapply(X, as.data.frame)
        stopifnot(`==`(lapply(X, nrow) %>% unlist %>% unique %>% length(), 1))
        
        res <- imap_dfr(X, ~{.getUpDown(.x) %>% mutate(block = .y)})
    }
    
    object <- list()
    object[["X"]] <- X
    object[["cluster"]] <- res
    class(object) <- "UpDown"
    return(object)
}

#' @importFrom plyr mapvalues
.getUpDown <- function(X){
    tmp <- lapply(X, function(x) {
        factor(sign(diff(x)), levels = c(1, -1, 0)) %>%
            plyr::mapvalues( from = c(1, -1, 0), to = c("Up", "Down", "0")) %>%
            as.character() %>%
            paste0(collapse = "_")}) %>% 
        as.data.frame() %>% t %>% as.data.frame() %>% rownames_to_column("molecule") %>%
        dplyr::rename("cluster"="V1")
    
    return(tmp)
}


# add getCluster for UpDown clusters
#' @export
getCluster.UpDown <- function(X, user.block = NULL, user.cluster = NULL){
    results <- X$cluster
    
    results <- filter.getCluster(X = results, user.block = user.block, user.cluster = user.cluster)
    class(results) <- c("cluster.df", "data.frame")
    return(results)
}
