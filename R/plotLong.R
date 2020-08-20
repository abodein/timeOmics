#' Plot Longitudinal Profiles by Cluster
#'
#' This function provides a expression profile representation over time and by cluster.
#'
#' @param object a mixOmics result of class (s)pca, (s)pls, block.(s)pls.
#' @param time (optional) a numeric vector, the same size as \code{ncol(X)}, to change the time scale.
#' @param plot a logical, if TRUE then  a plot is produced. Otherwise, the data.frame on which the plot is based on is returned.
#' @param center a logical value indicating whether the variables should be shifted to be zero centered. 
#' @param scale a logical value indicating whether the variables should be scaled to have unit variance before the analysis takes place. 
#' @param title character indicating the title plot.
#' @param X.label x axis titles.
#' @param Y.label y axis titles.
#' @param legend a logical, to display or not the legend.
#' @param legend.title if \code{legend} is provided, title of the legend.
#' @param legend.block.name a character vector corresponding to the size of the number of blocks in the mixOmics object. 
#' 
#'
#' @return
#' a data.frame (gathered form) containing the following columns:
#' \item{time}{x axis values}
#' \item{molecule}{names of features}
#' \item{value}{y axis values}
#' \item{cluster}{assigned clusters}
#' \item{block}{name of 'blocks'}
#' 
#' @seealso 
#' \code{\link[timeOmics]{getCluster}}
#'
#' @examples
#' demo <- suppressWarnings(get_demo_cluster())
#' X <- demo$X
#' Y <- demo$Y
#' Z <- demo$Z
#' 
#' # (s)pca
#' pca.res <- mixOmics::pca(X, ncomp = 3)
#' plotLong(pca.res)
#' spca.res <- mixOmics::spca(X, ncomp =2, keepX = c(15, 10))
#' plotLong(spca.res)
#' 
#' # (s)pls
#' pls.res <- mixOmics::pls(X,Y)
#' plotLong(pls.res)
#' spls.res <- mixOmics::spls(X,Y, keepX = c(15,10), keepY=c(5,6))
#' plotLong(spls.res)
#' 
#' # (s)block.spls
#' block.pls.res <- mixOmics::block.pls(X=list(X=X,Z=Z), Y=Y)
#' plotLong(block.pls.res)
#' block.spls.res <- mixOmics::block.spls(X=list(X=X,Z=Z), Y=Y, 
#'                              keepX = list(X = c(15,10), Z = c(5,6)), 
#'                              keepY = c(3,6))
#' plotLong(block.spls.res)
#' 
#'
#' @import ggplot2
#' @import mixOmics
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr mutate select left_join
#' @importFrom tidyr pivot_longer
#' @export
plotLong <- function(object, time = NULL, plot = TRUE, center = TRUE, scale = TRUE, 
                     title="Time-course Expression", X.label=NULL, Y.label=NULL, 
                     legend=FALSE, legend.title=NULL, legend.block.name = NULL)
{
    
    # Check parameters
    #-- object
    allowed_object = c("pca", "spca", "mixo_pls", "mixo_spls", "block.pls", "block.spls")
    if(!any(class(object) %in% allowed_object)){
        stop("invalid object, should be one of c(pca, spca, mixo_pls, mixo_spls, block.pls, block.spls)")
    }
    
    #-- plot : if plot is not correct, plot = FALSE
    if(is.null(plot)) plot = FALSE
    if(!is.finite(plot) || !is.logical(plot)){plot = FALSE}
    
    #-- center
    if(is.null(center)) center = TRUE
    if(!is.finite(center) || !is.logical(center)){center = TRUE}
    
    #-- scale
    if(is.null(scale)) scale = TRUE
    if(!is.finite(scale) || !is.logical(scale)){scale = TRUE}
    
    # graphical options
    #-- title
    if(!is.character(title)){title = NULL}
    
    #-- X.label
    if(!is.character(X.label)){X.label = NULL}
    
    #-- Y.label
    if(!is.character(Y.label)){Y.label = NULL}
    
    #-- legend
    if(is.null(legend)) legend = FALSE
    if(!is.finite(legend) || !is.logical(legend)){legend = FALSE}
    
    #-- legend title
    if(!is.character(legend.title)){legend.title = NULL}
    
    # cluster info
    cluster <- getCluster(object)
    
    #-- legend.block.name
    if(!is.null(legend.block.name)){
        check_legend.block.name(legend.block.name, cluster)
        new.block.name.tmp <- list(new.block = legend.block.name, block = unique(cluster$block)) %>%
            as.data.frame
        cluster <- cluster %>% left_join(new.block.name.tmp, by = c("block"="block")) %>%
            mutate(old.block = block, block = new.block) %>% 
            dplyr::select(-new.block)
    }
    
    if(is(object, "pca") || is(object, "spca")){
        #-- check time
        if(!is.null(time) && (!is_almostInteger_vector(time) || (length(time) != nrow(object$X)))){
            stop("'time' should be a numeric vector")
        }
        # scale/unscale if desired
        data <- unscale(object$X) %>%
            as.data.frame() %>%
            dplyr::select(intersect(cluster$molecule, colnames(.))) %>%
            scale(scale, center)
        
        
    } else if(is(object, "mixo_pls") || is(object, "mixo_spls")){
        #-- check time
        if(!is.null(time) && (!is_almostInteger_vector(time) || (length(time) != nrow(object$X)))){
            stop("'time' should be a numeric vector")
        }
        data.X <- unscale(object$X) %>%
            scale(scale, center)
        data.Y <- unscale(object$Y) %>%
            scale(scale, center)
        
        data <- cbind(data.X, data.Y) %>%
            as.data.frame() %>%
            dplyr::select(intersect(cluster$molecule, colnames(.)))
        
    } else if(is(object, "block.pls") || is(object, "block.spls")){
        #-- check time
        if(!is.null(time) &&(!is_almostInteger_vector(time) || (length(time) != nrow(object$X[[1]])))){
            stop("'time' should be a numeric vector")
        }
        
        data <- lapply(object$X, function(i){unscale(i) %>%
                scale(scale, center)}) %>% 
            do.call(what = "cbind")

        data <- as.data.frame(data) %>% 
            dplyr::select(intersect(cluster$molecule, colnames(.)))
    }
    
    # plot
    if(!is.null(time)){
        rownames(data) <- time   
    }
    data.gather <- data %>%
        as.data.frame() %>% 
        tibble::rownames_to_column("time") %>%
        dplyr::mutate(time = as.numeric(.$time)) %>%
        tidyr::pivot_longer(names_to = "molecule", values_to = "value", -time) %>%
        dplyr::left_join(cluster, by = c("molecule"="molecule")) %>%
        dplyr::mutate(block = factor(block))
    
    gg <- ggplot(data.gather, aes(x = time, y = value, group = molecule)) +
        geom_line(aes(color = block)) +
        facet_grid(contribution ~ comp, scales = "free") +
        scale_color_manual(values = mixOmics::color.mixo(1:length(levels(data.gather$block)))) +
        theme_bw() 
    
    if(is.character(title)){
        gg <- gg + ggtitle(title)
    }
    
    if(!is.null(X.label)){
        gg <- gg + xlab(X.label)
    } else {
        gg <- gg + xlab("Time")
    }
    
    if(!is.null(Y.label)){
        gg <- gg + ylab(Y.label)
    } else {
        gg <- gg + ylab("Expression")
    }
    
    if(!legend){
        gg <- gg + theme(legend.position = "none")
    } else { # legend is TRUE
        if(!is.null(legend.title)){
            gg <- gg + labs(color = legend.title)
        }
    }
    
    if(plot){
        print(gg)
    }
    return(invisible(gg$data))
}
