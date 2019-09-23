#' tuneCluster.spca
#'
#' @examples
#' demo <- suppressMessages(get_demo_cluster())
#' X <- demo$X
#' tune.spca.res <- tuneCluster.spca(X = X, ncomp = 2, test.keepX = c(2,5,7))
#' plot(tune.spca.res)
#' plot(tune.spca.res, comp = 2)

tuneCluster.spca <- function(X, ncomp = 2, test.keepX = rep(ncol(X), ncomp), ...){
    #-- checking input parameters ---------------------------------------------#
    #--------------------------------------------------------------------------#

    #-- X
    X <- check.matrix(X)
    if(!is.numeric.matrix(X)){
        stop("X must be a numeric matrix with finite value")
    }

    #-- ncomp
    if (is.null(ncomp))
        ncomp = min(nrow(X), ncol(X))

    if (!is.almostInteger(ncomp) || ncomp < 1)
        stop("invalid value for 'ncomp'.")
    if (ncomp > min(ncol(X), nrow(X)))
        stop("use smaller 'ncomp'")

    #-- keepX
    if (is.null(test.keepX) | length(test.keepX) == 1 | !is.numeric(test.keepX))
        stop("'test.keepX' must be a numeric vector with more than two entries")
    test.keepX <- sort(unique(test.keepX))

    min.test.keepX <- rep(min(test.keepX), ncomp)


    #-- launch tuning  --------------------------------------------------------#
    #--------------------------------------------------------------------------#

    #-- 0. set output object
    result = list()

    #-- 1. compute dissimilarity matrix for silhouette coef. (once and for all)
    dmatrix <- dmatrix.spearman.dissimilarity(X)
    cluster <- as.data.frame(list("feature" = rownames(dmatrix)))

    for(comp in 1:ncomp){
        result[[comp]] <- list()

        tmp.keepX <- min.test.keepX  # foreach comp, keepX of other comp is set to minimum
        for(keepX in test.keepX){
            result[[comp]][[as.character(keepX)]] <- list()
            tmp.keepX[comp] <- keepX

            #-- 2. run spca
            kX = tmp.keepX
            spca.res <- mixOmics::spca(X = X, ncomp = ncomp, keepX = kX)

            #-- 3. extract clusters
            tmp.cluster <- getCluster(spca.res)
            tmp.cluster <- suppressWarnings(dplyr::left_join(cluster, tmp.cluster[c(1,4)],
                                                             by = c("feature"="molecule"))) %>%
                mutate(cluster = as.numeric(as.character(cluster))) %>%
                mutate(cluster = ifelse(is.na(cluster), 0, cluster))

            #-- 4. compute silhouette
            sil <- silhouette(dmatrix, tmp.cluster$cluster)
            result[[comp]][[as.character(keepX)]]$average.cluster <- sil$average.cluster
            result[[comp]][[as.character(keepX)]]$average <- sil$average
            result[[comp]][[as.character(keepX)]]$average.n0 <- sil$feature %>%
                filter(cluster != 0) %>% pull(silhouette.coef) %>% mean
        }
    }
    class(result) <- "spca.tune.silhouette"
    return(result)
}

plot.spca.tune.silhouette <- function(X, comp = 1, plot = TRUE){
    #-- checking input parameters ---------------------------------------------#
    #--------------------------------------------------------------------------#

    #-- should be a spca.tune.silhouette" object.

    #-- comp
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
