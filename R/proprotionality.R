#' Proportionality Distance
#' 
#' \code{proportionality} is a wrapper that compute proportionality distance for 
#' a clustering result (\code{pca}, \code{spca}, \code{pls}, \code{spls}, \code{block.pls}, \code{block.spls}).
#' and it performs a u-test to compare the median within a cluster to the median of the entire background set.
#' 
#' @param X an object of the class: \code{pca}, \code{spca}, \code{pls}, \code{spls}, \code{block.pls} or \code{block.spls}
#' 
#' @return 
#' Return a list containing the following components:
#'   \item{propr.distance}{Square matrix with proportionality distance between pairs of features}
#'   \item{propr.distance.w.cluster}{distance between pairs with cluster label}
#'   \item{pvalue}{Wilcoxon U-test p-value comparing the medians within clusters and with the entire background set}
#' 
#' @seealso 
#' \code{\link[propr]{propr}}
#' 
#' @references 
#' Lovell, D., Pawlowsky-Glahn, V., Egozcue, J. J., Marguerat, S., Bähler, J. (2015). Proportionality: a valid alternative to correlation for relative data. PLoS Comput. Biol. 11, e1004075. doi: 10.1371/journal.pcbi.1004075
#' 
#' Quinn, T. P., Richardson, M. F., Lovell, D., Crowley, T. M. (2017). propr: an r-package for identifying proportionally abundant features using compositional data analysis. Sci. Rep. 7, 16252. doi: 10.1038/s41598-017-16520-0
#' 
#' @examples 
#' demo <- get_demo_cluster()
#' 
#' # pca
#' X <- demo$pca
#' propr.res <- proportionality(X)
#' plot(propr.res)
#' 
#' # pls
#' X <- demo$spls
#' propr.res <- proportionality(X)
#' plot(propr.res)
#' 
#' # block.pls
#' X <- demo$block.spls
#' propr.res <- proportionality(X)
#' plot(propr.res)
#' 
#' @importFrom dplyr mutate filter rename left_join
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @importFrom propr propr
#' 
#' @export
proportionality <- function(X){
    stopifnot(is(X, c("pca", "spca", "mixo_pls", "mixo_spls", "block.pls", "block.spls")))
    
    # 1. get cluster
    cluster.info <- getCluster(X) %>% dplyr::select(molecule, cluster)
    
    # 2. extract data and add cluster
    
    # pca / spca
    if(is(X,c("pca", "spca"))){
        # unscaling + positive value
        data <- unscale(X$X) %>% `+`(abs(min(.))) 
    } else if(is(X, c("mixo_pls", "mixo_spls"))){
        # unscale X, unscale Y, cat
        data.X <- unscale(X$X) %>% `+`(abs(min(.)))
        data.Y <- unscale(X$Y) %>% `+`(abs(min(.)))
        data <- cbind(data.X, data.Y)
    } else if(is(X, c("block.pls", "block.spls"))){
        if(is.null(X$Y)){
            data <- lapply(X$X, function(x) x %>% unscale %>% `+`(abs(min(.)))) %>%
                do.call(what="cbind")
        } else {
            data.X <- lapply(X$X, function(x) x %>% unscale %>% `+`(abs(min(.)))) %>%
                do.call(what="cbind")
            data.Y <- unscale(X$Y) %>% `+`(abs(min(.)))
            data <- cbind(data.X, data.Y)
        }
    }
    
    # 3. compute phi_s
    data.propr <- suppressMessages(as.data.frame(propr::propr(data, metric = "phs")@matrix))
    # gather, add cluster info, compute basic stats
    data.propr.gather <- data.propr %>% tibble::rownames_to_column("feature1") %>%
        tidyr::pivot_longer(-feature1, names_to = "feature2", values_to = 'value') %>%
        dplyr::filter(feature1 %in% cluster.info$molecule) %>%
        dplyr::left_join(cluster.info, by = c("feature1"="molecule")) %>%
        dplyr::rename("cluster1"="cluster") %>%
        dplyr::filter(feature2 %in% cluster.info$molecule) %>%
        dplyr::left_join(cluster.info, by = c("feature2"="molecule")) %>%
        dplyr::rename("cluster2"="cluster") %>%
        dplyr::mutate(insideout = ifelse(cluster1 == cluster2, "inside", "outside"))
        
    # compute stat (median, u-test pval, adj.pval)
    res.stat <- stat_median(data.propr.gather) %>% na.omit()
    
    res <- list()
    res[["propr.distance"]] <- data.propr
    res[["propr.distance.w.cluster"]] <- data.propr.gather
    res[["pvalue"]] <- res.stat
    class(res) <- "proportionality"
    return(res)
}

#' @importFrom dplyr mutate filter pull
#' @importFrom purrr set_names
#' @importFrom magrittr %>%
stat_median <- function(res.phs.X){
    i = 1
    res.pval <- matrix(ncol = 4, nrow = 4) %>% as.data.frame() %>%
        purrr::set_names("cluster", "median_inside", "median_outside", "Pvalue")
    for(clu in unique(res.phs.X$cluster1)){
        inside <- res.phs.X %>% filter(cluster1 == clu) %>% 
            dplyr::filter(cluster2==clu) %>% 
            dplyr::pull(value)
        outside <- res.phs.X %>% 
            filter(cluster1 == clu) %>% 
            dplyr::filter(cluster2!=clu) %>% 
            dplyr::pull(value)
        
        #ttest.pval <- t.test(inside, outside)$p.value
        utest.pval <- stats::wilcox.test(inside, outside)$p.value
        
        res.pval[i,] <- c(clu, round(median(inside), digits = 2),
                          round(median(outside), digits = 2), utest.pval)
        i = i+1
    }
    as.data.frame(res.pval) %>% 
        dplyr::mutate("Adj.Pvalue" = stats::p.adjust(Pvalue, method = "fdr")) %>%
        na.omit
    return(res.pval)
}

#' @export
#' @import ggplot2
#' @importFrom mixOmics color.mixo
plot.proportionality <- function(x, ...){
    ggplot2::ggplot(data = x$propr.distance.w.cluster, 
                    aes(x=as.factor(cluster1), y=value, col=insideout)) + 
        geom_boxplot() + theme_bw() + 
        xlab("Cluster ID") + 
        ylab("Proportionality distance") + 
        labs(color = "Proportionality distance") +
        scale_color_manual(values = mixOmics::color.mixo(1:2))
}