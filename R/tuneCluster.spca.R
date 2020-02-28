#' tuneCluster.spca
#'
#' This function identify the number of feautures to keep per component and thus by cluster in \code{mixOmics::spca} 
#' by optimizing the silhouette coefficient, which assesses the quality of clustering.
#' 
#' @param X numeric matrix (or data.frame) with features in columns and samples in rows
#' @param ncomp integer, number of component to include in the model
#' @param test.keepX vector of integer containing the different value of keepX to test for block \code{X}.
#' @param ... other parameters to be included in the spls model (see \code{mixOmics::spca})
#' 
#' @return 
#' \item{silhouette}{silhouette coef. computed for every combinasion of keepX/keepY}
#' \item{ncomp}{number of component included in the model}
#' \item{test.keepX}{list of tested keepX}
#' \item{block}{names of blocks}
#' \item{slopes}{"slopes" computed from the silhouette coef. for each keepX and keepY, used to determine the best keepX and keepY}
#' \item{choice.keepX}{best \code{keepX} for each component}
#' 
#' @details
#' For each component and for each keepX value, a spls is done from these parameters.
#' Then the clustering is performed and the silhouette coefficient is calculated for this clustering.
#'
#' We then calculate "slopes" where keepX are the coordinates and the silhouette is the intensity.
#' A z-score is assigned to each slope.
#' We then identify the most significant slope which indicates a drop in the silhouette coefficient and thus a deterioration of the clustering.
#'
#' 
#' @examples
#' demo <- get_demo_cluster()
#' X <- demo$X
#' 
#' # tuning
#' tune.spca.res <- tuneCluster.spca(X = X, ncomp = 2, test.keepX = c(2:10))
#' keepX <- tune.spca.res$choice.keepX
#' 
#' # final model
#' spca.res <- spca(X=X, ncomp = 2, keepX = keepX)
#' plot(spca.res)
#' 


#' @import mixOmics
#' @importFrom dplyr left_join
#' @importFrom dplyr mutate
#' @importFrom dplyr filter
#'
#' @export
tuneCluster.spca <- function(X, ncomp = 2, test.keepX = rep(ncol(X), ncomp), ...){
    #-- checking input parameters ---------------------------------------------#
    #--------------------------------------------------------------------------#

    #-- X
    X <- validate.matrix.X(X)

    #-- ncomp
    ncomp <- validate.ncomp(ncomp, list(X))
    
    #-- test.keepX
    test.keepX <- validate.test.keepX(test.keepX = test.keepX, X = X)
    min.test.keepX <- rep(min(test.keepX), ncomp)


    #-- launch tuning  --------------------------------------------------------#
    #--------------------------------------------------------------------------#

    #-- 0. set output object
    result <- as.data.frame(matrix(ncol = 4, nrow = length(test.keepX)*ncomp))
    colnames(result) <- c("comp", "X", "pos", "neg")
    result.index <- 0

    #-- 1. compute dissimilarity matrix for silhouette coef. (once and for all)
    dmatrix <- dmatrix.spearman.dissimilarity(X)
    cluster <- as.data.frame(list("feature" = rownames(dmatrix)))

    #-- tuning
    for(comp in 1:ncomp){
        tmp.keepX <- min.test.keepX  # foreach comp, keepX of other comp is set to minimum
        for(keepX in test.keepX){
            tmp.keepX[comp] <- keepX

            #-- 2. run spca
            kX = tmp.keepX
            spca.res <- mixOmics::spca(X = X, ncomp = ncomp, keepX = kX)

            #-- 3. extract clusters
            tmp.cluster <- getCluster(spca.res)
            tmp.cluster <- suppressWarnings(dplyr::left_join(cluster, tmp.cluster[c(1,4)],
                                                             by = c("feature"="molecule"))) %>%
                dplyr::mutate(cluster = as.numeric(as.character(cluster))) %>%
                dplyr::mutate(cluster = ifelse(is.na(cluster), 0, cluster))

            #-- 4. compute silhouette
            sil <- silhouette(dmatrix, tmp.cluster$cluster)
            
            #-- 6. store
            result.index <- result.index + 1
            result[result.index, "comp"] <- comp
            result[result.index, "X"] <- kX[comp]
            
            pos.res <-  sil$average.cluster  %>% 
                dplyr::filter(cluster == comp) %>% dplyr::pull(silhouette.coef)
            result[result.index, "pos"] <- ifelse(length(pos.res) == 0, NA, pos.res)
            neg.res <-  sil$average.cluster  %>%
                dplyr::filter(cluster == -comp) %>% dplyr::pull(silhouette.coef)
            result[result.index, "neg"] <- ifelse(length(neg.res) == 0, NA, neg.res)
        }
    }
    result <- list("silhouette" = result)
    result[["ncomp"]] <- ncomp
    result[["test.keepX"]] <- test.keepX
    result[["block"]] <- c("X")
    class(result) <- "spca.tune.silhouette"

    #-- 7. choice.keepX
    result[["slopes"]] <- tune.silhouette.get_slopes(result)
    tmp <- tune.silhouette.get_choice_keepX(result) # choice keepX/keepY
    result[["choice.keepX"]] <- unlist(lapply(tmp, function(x) x$X))
    return(result)
}

#' plot.spca.tune.silhouette
#'
#' Plot 
#'
#' @import mixOmics
#' @importFrom dplyr left_join mutate filter
#' @import ggplot2
#' @export
plot.spca.tune.silhouette <- function(X, comp = 1, plot = TRUE){
    #-- checking input parameters ---------------------------------------------#
    #--------------------------------------------------------------------------#

    #-- should be a spca.tune.silhouette" object.

    #-- comp
    if(length(comp) != 1){
        stop(paste0("Invalid 'comp', shoud be an integer between 1 and ", length(X)))
    }
    if(!(comp %in% seq_along(X))){
        stop(paste0("Invalid 'comp', shoud be an integer between 1 and ", length(X)))
    }
    

    #-- plot : if plot is not correct, plot = FALSE
    if(is.null(plot)) plot = FALSE
    if(!is.finite(plot) || !is.logical(plot)){ plot = FALSE }

    ncomp <- length(X)
    test.keepX <- names(X[[1]])
    X.df <- .spca.tune.rearrange_result(X)

    plot.df <- X.df[[comp]] %>% filter(cluster != 0) %>%
        mutate(keepX = as.numeric(keepX)) %>%
        mutate(cluster = as.numeric(as.character(cluster))) %>%
        mutate(linetype = factor(ifelse(sign(cluster) > 0, "pos", "neg"), levels = c("pos", "neg"))) %>%
        mutate(alpha = I(ifelse(abs(cluster) == comp, 1, 0.5))) %>%
        mutate(size = I(ifelse(abs(cluster) == comp, 1, 0.5))) %>%
        mutate(color = paste("comp", abs(cluster)))

    ggplot.df <- ggplot(plot.df, aes(x = keepX, y = silhouette.coef, color = color, group = cluster)) +
        geom_line(aes(alpha = alpha, size = size, lty = linetype)) +
        scale_x_continuous(name = "keepX", breaks = as.numeric(test.keepX), labels = test.keepX)+
        scale_color_manual(values = color.mixo(1:ncomp))  +
        facet_wrap(~paste0("comp ", comp)) +
        theme_bw() + labs(y = "Silhouette Coefficient", color = "Comp.", lty = "Contrib.")

    if(plot){
        print(ggplot.df)
    }
    return(invisible(ggplot.df))
}

.spca.tune.rearrange_result <- function(X){
    res <- lapply(seq_along(X), function(Z){
        do.call("rbind",imap(X[[Z]], ~.x$average.cluster %>%
                                 mutate(keepX = .y))) %>%
            mutate(comp = Z)})
}
