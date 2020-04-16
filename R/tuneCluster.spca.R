#' Feature Selection Optimization for sPCA method 
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
#' demo <- suppressWarnings(get_demo_cluster())
#' X <- demo$X
#' 
#' # tuning
#' tune.spca.res <- tuneCluster.spca(X = X, ncomp = 2, test.keepX = c(2:10))
#' keepX <- tune.spca.res$choice.keepX
#' plot(tune.spca.res)
#' 
#' # final model
#' spca.res <- mixOmics::spca(X=X, ncomp = 2, keepX = keepX)
#' plotLong(spca.res)


#' @import mixOmics
#' @importFrom dplyr left_join
#' @importFrom dplyr mutate
#' @importFrom dplyr filter
#' @export
tuneCluster.spca <- function(X, ncomp = 2, test.keepX = rep(ncol(X), ncomp), ...){
    #-- checking input parameters ---------------------------------------------#
    #--------------------------------------------------------------------------#

    #-- X
    X <- validate_matrix_X(X)

    #-- ncomp
    ncomp <- validate_ncomp(ncomp, list(X))
    
    #-- test.keepX
    test.keepX <- validate_test_keepX(test.keepX = test.keepX, X = X)
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
                dplyr::filter(cluster == comp) %>%
                dplyr::pull(silhouette.coef)
            result[result.index, "pos"] <- ifelse(length(pos.res) == 0, NA, pos.res)
            neg.res <-  sil$average.cluster  %>%
                dplyr::filter(cluster == -comp) %>%
                dplyr::pull(silhouette.coef)
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




#' @import mixOmics
#' @importFrom dplyr left_join mutate filter group_by top_n
#' @import ggplot2
#' @importFrom tidyr pivot_longer
#' @importFrom ggrepel geom_label_repel
#' @export
plot.spca.tune.silhouette <- function(x, ...){
    #-- checking input parameters ---------------------------------------------#
    #--------------------------------------------------------------------------#

    #-- should be a spca.tune.silhouette" object.
    ncomp <- x$ncomp
    test.keepX <- x$test.keepX
    
    tmp <- x$silhouette %>% 
        tidyr::pivot_longer(names_to = "contrib", values_to = "value", -c(comp,X)) %>%
        na.omit %>% # remove NA intruced in pos/neg
        dplyr::mutate(comp = as.factor(paste0("Comp ", comp)), 
               contrib = factor(contrib, levels = c("pos","neg"))) %>%
        dplyr::group_by(comp, contrib)
    
    choice <- list(comp= as.factor(paste0("Comp ",names(x$choice.keepX))), 
                   X = unname(x$choice.keepX)) %>%
        as.data.frame() %>%
        dplyr::left_join(tmp, by = c("comp"="comp", "X"="X")) %>% 
        dplyr::group_by(comp, X) %>%
        dplyr::top_n(n = 1, wt = value)
    
    choice.vline <- choice %>%
        dplyr::select(c("comp", "X"))
    
    ggplot.df <- ggplot(tmp, aes(x=X, y =value, col = comp)) + 
        geom_line(aes(lty = contrib)) + facet_wrap(~comp) +
        theme_bw() +
        geom_point(data = choice, size = 5, pch = 18) +
        ggrepel::geom_label_repel(data = choice, aes(label = X), col = "black") +
        scale_color_manual(values = mixOmics::color.mixo(1:x$ncomp)) +
        labs(x ="tested keepX",  y = "Silhouette Coef.", color = "Comp.", lty = "Contrib.") +
        geom_vline(data = choice.vline, aes(xintercept = X), lty = 5, col = "grey")


    print(ggplot.df)
    
    return(invisible(ggplot.df))
}