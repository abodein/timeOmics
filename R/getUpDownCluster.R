#' Up-Down clustering
#' 
#' Performs a clustering based on the signs of variation between 2 timepoints.
#' Optionally, if the difference between 2 timepoints is lower than a given threshold, 
#' the returned difference will be 0.
#' 
#' @param X a dataframe or list of dataframe with the same number of rows.
#' @param diff_threshold a number (optional, default 0), if the difference between 2 values is lower than the threshold, the returned sign will be 0 (no variation).
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
#' res <- getUpDownCluster(X, diff_thresold = 15)
#' res_cluster <- getCluster(res)

#' @importFrom purrr imap_dfr
#' @importFrom checkmate check_number
#' 
#' @export
getUpDownCluster <- function(X, diff_thresold = 0){

    #stopifnot(class(X) %in% c("matrix", "data.frame", "list"))
    stopifnot(is(X, "matrix") || is(X, "data.frame") || is(X, "list"))
    checkmate::check_number(diff_thresold, null.ok = TRUE)
    
    
    if(is.matrix(X) || is.data.frame(X)){
        
        # check X
        X <- validate_matrix_X(X)
        X <- as.data.frame(X)
        
        res <- .getUpDown(X, diff_thresold = diff_thresold) %>% mutate(block = "X")
    }
    else if(is.list(X) & length(X)>1){
        
        # check X list
        X <- validate_list_matrix_X(X)
        X <- lapply(X, as.data.frame)
        stopifnot(`==`(lapply(X, nrow) %>% unlist %>% unique %>% length(), 1))
        
        res <- imap_dfr(X, ~{.getUpDown(.x, diff_thresold = diff_thresold) %>% mutate(block = .y)})
    }
    
    object <- list()
    object[["X"]] <- X
    object[["cluster"]] <- res
    class(object) <- "UpDown"
    return(object)
}

#' @importFrom plyr mapvalues
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr rename
.getUpDown <- function(X, diff_thresold){
    tmp <- lapply(X, function(x) {
        factor(sign(.apply_fc_threshold(diff(x), diff_thresold = diff_thresold)),
               levels = c(1, -1, 0)) %>%
            plyr::mapvalues( from = c(1, -1, 0), to = c("Up", "Down", "0")) %>%
            as.character() %>%
            paste0(collapse = "_")})
    tmp <- as.data.frame(tmp, check.names = FALSE) %>% 
        t %>% as.data.frame(check.names = FALSE) %>% 
        tibble::rownames_to_column("molecule") %>%
        dplyr::rename("cluster"="V1")
    return(tmp)
}


#' @examples 
#' demo <- suppressWarnings(get_demo_cluster())
#' x <- diff(demo$X[,1])
#' diff_thresold <- 15
.apply_fc_threshold <- function(x, diff_thresold){
    # x is numeric from diff function
    # threshold is numeric
    res <-  ifelse(abs(x) < diff_thresold, 0, x)
    return(res)
}

# add getCluster for UpDown clusters
#' @export
getCluster.UpDown <- function(X, user.block = NULL, user.cluster = NULL){
    results <- X$cluster
    
    results <- filter.getCluster(X = results, user.block = user.block, user.cluster = user.cluster)
    class(results) <- c("cluster.df", "data.frame")
    return(results)
}
