#' plotLong
#'
#' A more detailed description.
#'
#' @param X A numeric matrix.
#'
#' @param (s)pca, (s)pls, block.(s)pls results
#' that contains clustering information with molecule and cluster
#'
#'
#' @return
#' \describe{
#'   \item{One}{First item}
#'   \item{Two}{Second item}
#' }
#'
#' @examples
#'
#'
#' @export
plotLong <- function(object, time = NULL, plot = TRUE, center = TRUE, scale = TRUE){
    #allowed_object = c("pca", "spca", "mixo_pls", "mixo_spls", "block.pls", "block.spls")
    #stopifnot(is(object, allowed_object))
    UseMethod("plotLong")
}

#' @importFrom magrittr %>%
#' @export
plotLong.pca <- function(object, time = NULL, plot = TRUE, center = TRUE, scale = TRUE){

    if(!is.null(time)){
        stopifnot(is.numeric(time)) # works with integer
        stopifnot(length(time) == nrow(object$X))

    } else{ # time IS NULL
        # we rely on rownames and assume it's numeric,  correspond to times
        time <- as.numeric(rownames(object$X))
        # not numeric value can introduce NA
        stopifnot(!is.na(time))
    }
    # unscale and rescale if desired
    data <- object$X %>% scale(scale, center)

    # cluster info
    cluster <- getCluster(object)

    gg <- plotLongGGplot(data = data, time = time, cluster = cluster)
    return(invisible(gg))
}

#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @export
plotLong.spca <- function(object, time = NULL, plot = TRUE, center = TRUE, scale = TRUE){
    print("plotLong.spca")

    if(!is.null(time)){
        stopifnot(is.numeric(time)) # works with integer
        stopifnot(length(time) == nrow(object$X))

    } else{ # time IS NULL
        # we rely on rownames and assume it's numeric,  correspond to times
        time <- as.numeric(rownames(object$X))
        # not numeric value can introduce NA
        stopifnot(!is.na(time))
    }

    # cluster info
    # unselected features are removed here.
    cluster <- getCluster(object)

    # unscale and rescale if desired
    data <- unscale(object$X) %>% as.data.frame() %>%
        dplyr::select(cluster$molecule) %>%
        scale(scale, center)


    gg <- plotLongGGplot(data = data, time = time, cluster = cluster)
    return(invisible(gg))
}

#' @export
#' @importFrom magrittr %>%
plotLong.mixo_pls <- function(object, time = NULL, plot = TRUE, center = TRUE, scale = TRUE){
    print("plotLong.mixo_pls")

    stopifnot(!is.null(object$Y)) # wrong object
    stopifnot(nrow(object$X) == nrow(object$Y)) # PLS so X and Y should
    # have same nrow.
    stopifnot(rownames(object$X) == rownames(object$Y))

    if(!is.null(time)){
        stopifnot(is.numeric(time)) # works with integer
        stopifnot(length(time) == nrow(object$X)) # can use object$Y

    } else{ # time IS NULL
        # we rely on rownames and assume it's numeric,  correspond to times
        time <- as.numeric(rownames(object$X))  # idem can use object$Y
        # not numeric value can introduce NA
        stopifnot(!is.na(time))
    }
    # unscale and rescale if desired
    # X
    X <- unscale(object$X) %>% scale(scale, center)
    Y <- unscale(object$Y) %>% scale(scale, center)
    data <- cbind(X,Y)

    # cluster info
    cluster <- getCluster(object)

    gg <- plotLongGGplot(data = data, time = time, cluster = cluster)
    return(invisible(gg))
}

#' @export
#' @importFrom magrittr %>%
#' @importFrom dplyr select
plotLong.mixo_spls <- function(object, time = NULL, plot = TRUE, center = TRUE, scale = TRUE){
    print("plotLong.mixo_spls")

    stopifnot(!is.null(object$Y)) # wrong object
    stopifnot(nrow(object$X) == nrow(object$Y)) # PLS so X and Y should
    # have same nrow.
    stopifnot(rownames(object$X) == rownames(object$Y))

    if(!is.null(time)){
        stopifnot(is.numeric(time)) # works with integer
        stopifnot(length(time) == nrow(object$X)) # can use object$Y

    } else{ # time IS NULL
        # we rely on rownames and assume it's numeric,  correspond to times
        time <- as.numeric(rownames(object$X))  # idem can use object$Y
        # not numeric value can introduce NA
        stopifnot(!is.na(time))
    }
    # unscale and rescale if desired
    # X
    X <- unscale(object$X) %>% scale(scale, center)
    Y <- unscale(object$Y) %>% scale(scale, center)

    # cluster info
    # unselected features are removed here.
    cluster <- getCluster(object)


    data <- cbind(X,Y) %>% as.data.frame %>%
        dplyr::select(cluster$molecule)

    gg <- plotLongGGplot(data = data, time = time, cluster = cluster)
    return(invisible(gg))
}

#' @export
#' @importFrom magrittr %>%
plotLong.block.pls <- function(object, time = NULL, plot = TRUE, center = TRUE, scale = TRUE){
    print("plotLong.block.pls")

    stopifnot(1 == object$X %>% lapply(nrow) %>% unlist %>% unique %>% length())
    # have same nrow.
    stopifnot(lapply(object$X, function(x)
        rownames(x) == rownames(object$X[[1]])) %>%
            unlist %>% all) # have same rownames

    if(!is.null(time)){
        stopifnot(is.numeric(time)) # works with integer
        stopifnot(length(time) == nrow(object$X[[1]])) # can use object$Y

    } else{ # time IS NULL
        # we rely on rownames and assume it's numeric,  correspond to times
        time <- as.numeric(rownames(object$X[[1]]))  # idem can use object$Y
        # not numeric value can introduce NA
        stopifnot(!is.na(time))
    }
    # unscale and rescale if desired
    data <- lapply(object$X, function(x) unscale(x) %>% scale(scale, center))

    # cluster info
    # unselected features are removed here.
    cluster <- getCluster(object)


    data <- do.call("cbind", data) %>% as.data.frame

    gg <- plotLongGGplot(data = data, time = time, cluster = cluster)
    return(invisible(gg))
}

#' @export
#' @importFrom magrittr %>%
#' @importFrom dplyr select
plotLong.block.spls <- function(object, time = NULL, plot = TRUE, center = TRUE, scale = TRUE){
    print("plotLong.block.spls")

    stopifnot(1 == object$X %>% lapply(nrow) %>% unlist %>% unique %>% length())
    # have same nrow.
    stopifnot(lapply(object$X, function(x)
        rownames(x) == rownames(object$X[[1]])) %>%
            unlist %>% all) # have same rownames

    if(!is.null(time)){
        stopifnot(is.numeric(time)) # works with integer
        stopifnot(length(time) == nrow(object$X[[1]])) # can use object$Y

    } else{ # time IS NULL
        # we rely on rownames and assume it's numeric,  correspond to times
        time <- as.numeric(rownames(object$X[[1]]))  # idem can use object$Y
        # not numeric value can introduce NA
        stopifnot(!is.na(time))
    }
    # unscale and rescale if desired
    data <- lapply(object$X, function(x) unscale(x) %>% scale(scale, center))

    # cluster info
    # unselected features are removed here.
    cluster <- getCluster(object)


    data <- do.call("cbind", data) %>% as.data.frame %>%
        dplyr::select(cluster$molecule)

    gg <- plotLongGGplot(data = data, time = time, cluster = cluster)
    return(invisible(gg))
}

#' @export
#' @importFrom mixOmics color.mixo
#' @import ggplot2
#' @importFrom tidyr gather
#' @importFrom dplyr mutate left_join
plotLongGGplot <- function(data, time, cluster, plot = TRUE){
    # graphical call
    stopifnot(is.numeric(time)) # works with integer
    stopifnot(nrow(data) == length(time))

    data.gather <- data %>% as.data.frame() %>%
        mutate(time = time) %>%
        gather(molecule, value, -time) %>%
        left_join(cluster, by = c("molecule"="molecule"))

    gg <- ggplot(data.gather, aes(x = time, y = value, group = molecule)) +
        geom_line(aes(color = block)) +
        facet_grid(contribution ~ comp, scales = "free") +
        scale_color_manual(values = mixOmics::color.mixo(1:length(levels(data.gather$cluster)))) +
        theme_bw()

    if(plot){
        print(gg)
    }
    return(data.gather)
}
