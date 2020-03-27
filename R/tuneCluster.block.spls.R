#' Feature Selection Optimization for block (s)PLS method 
#'
#' This function identify the number of feautures to keep per component and thus by cluster in \code{mixOmics::block.spls} by optimizing the silhouette coefficient, which assesses the quality of clustering.
#' 
#' @param X list of numeric matrix (or data.frame) with features in columns and samples in rows (with samples order matching in all data sets).
#' @param Y (optional) numeric matrix (or data.frame) with features in columns and samples in rows (same rows as \code{X}).
#' @param indY integer, to supply if Y is missing, indicates the position of the matrix response in the list \code{X}.
#' @param ncomp integer, number of component to include in the model
#' @param test.list.keepX list of integers with the same size as X. Each entry corresponds to the different keepX value to test for each block of \code{X}.
#' @param test.keepY only if Y is provideid. Vector of integer containing the different value of keepY to test for block \code{Y}.
#' @param ... other parameters to be included in the spls model (see \code{mixOmics::block.spls})
#' 
#' @return 
#' \item{silhouette}{silhouette coef. computed for every combinasion of keepX/keepY}
#' \item{ncomp}{number of component included in the model}
#' \item{test.keepX}{list of tested keepX}
#' \item{test.keepY}{list of tested keepY}
#' \item{block}{names of blocks}
#' \item{slopes}{"slopes" computed from the silhouette coef. for each keepX and keepY, used to determine the best keepX and keepY}
#' \item{choice.keepX}{best \code{keepX} for each component}
#' \item{choice.keepY}{best \code{keepY} for each component}
#' 
#'
#' @details
#' For each component and for each keepX/keepY value, a spls is done from these parameters.
#' Then the clustering is performed and the silhouette coefficient is calculated for this clustering.
#'
#' We then calculate "slopes" where keepX/keepY are the coordinates and the silhouette is the intensity.
#' A z-score is assigned to each slope.
#' We then identify the most significant slope which indicates a drop in the silhouette coefficient and thus a deterioration of the clustering.
#'
#' 
#' @seealso 
#' \code{\link[mixOmics]{block.spls}}, \code{\link[timeOmics]{getCluster}}, \code{\link[timeOmics]{plotLong}}
#'
#' @examples
#' demo <- get_demo_cluster()
#' X <- list(X = demo$X, Z = demo$Z)
#' Y <- demo$Y
#' test.list.keepX <- list("X" = c(5,10,15,20), "Z" = c(2,4,6,8))
#' test.keepY <- c(2:5)
#' 
#' # tuning
#' tune.block.spls <- tuneCluster.block.spls(X= X, Y= Y, 
#'                                           test.list.keepX= test.list.keepX, 
#'                                           test.keepY= test.keepY, 
#'                                           mode= "canonical")
#' keepX <- tune.block.spls$choice.keepX
#' keepY <- tune.block.spls$choice.keepY
#' 
#' # final model
#' block.spls.res <- mixOmics::block.spls(X= X, Y= Y, keepX = keepX, 
#'                              keepY = keepY, ncomp = 2, mode = "canonical")
#' # get clusters and plot longitudinal profile by cluster
#' block.spls.cluster <- getCluster(block.spls.res)

#' @import mixOmics
#' @importFrom dplyr left_join
#' @importFrom dplyr mutate
#' @importFrom dplyr filter
#' @export
tuneCluster.block.spls <- function(X, Y = NULL,  indY = NULL, ncomp = 2, 
                                   test.list.keepX = NULL, test.keepY = NULL, ...){
    #-- checking input parameters ---------------------------------------------#
    #--------------------------------------------------------------------------#

    #-- X
    X <- validate.list.matrix.X(X)

    #-- Y
    if(!is.null(Y)){
        Y <- validate.matrix.Y(Y)
    } else {
        indY <- validate.indY(indY, X)
    }

    #-- ncomp
    ncomp <- validate.ncomp(ncomp = ncomp, X =X)

    #-- keepX
    test.list.keepX <- validate.test.list.keepX(test.keepX = test.list.keepX, X = X, ncomp = ncomp)

    #-- keepY
    test.keepY <- validate.test.keepY(test.keepY = test.keepY, Y = Y)

    list.keepX.keepY <- test.list.keepX
    if(!is.null(Y)) {
        list.keepX.keepY$Y <- test.keepY
    }
    list.keepX.keepY <- expand.grid(list.keepX.keepY, stringsAsFactors = FALSE,
                                    KEEP.OUT.ATTRS = FALSE)


    #-- launch tuning  --------------------------------------------------------#
    #--------------------------------------------------------------------------#

    #-- 0. set output object
    result <- matrix(ncol = 3 + length(X),nrow = nrow(list.keepX.keepY)*ncomp) %>%
        as.data.frame() %>% purrr::set_names(c("comp", names(X), "pos", "neg"))
    if(!is.null(Y)){
        result <- result %>% mutate("Y" = NA)
    }
    result.index <- 1

    #--1. compute dissimilarity matrix for silhouette coef. (once and for all)
    all_data <- X
    all_data_name <- as.character(unlist(imap(X, ~rep(.y,ncol(.x)))))
    if(is.null(all_data$Y) & !is.null(Y)){
        all_data$Y <- Y
        all_data_name <- c(all_data_name, rep("Y", ncol(Y)))
    }
    all_data <- do.call("cbind", all_data)
    dmatrix <- dmatrix.spearman.dissimilarity(all_data)

    cluster <- as.data.frame(list("feature" = rownames(dmatrix),
                                  "block" = all_data_name))

    #--2. tuning
    for(comp in 1:ncomp){
        for(index.list.kX.kY in 1:nrow(list.keepX.keepY)){
            # foreach comp, keepX and keepY of other comp is set to minimum
            kX <- lapply(list.keepX.keepY, function(x) rep(min(x), ncomp))
            for(block in names(list.keepX.keepY)){
                kX[[block]][comp] <- list.keepX.keepY[index.list.kX.kY, block]
                result[result.index, block] = list.keepX.keepY[index.list.kX.kY, block]
            }
            if(!is.null(Y)){
                kY <- kX[["Y"]]
                kX[["Y"]] <- NULL

                #--3. run spls
                block.spls.res <- mixOmics::block.spls(X = X, Y = Y, ncomp =  ncomp,
                                                       keepX = kX, keepY = kY, ...)
            } else {
                #--3. run spls
                block.spls.res <- mixOmics::block.spls(X = X, ncomp =  ncomp,
                                                       keepX = kX, indY = indY, ...)
            }


            #--4. extract clusters
            tmp.cluster <- getCluster(block.spls.res)
            tmp.cluster <- suppressWarnings(dplyr::left_join(cluster, tmp.cluster[c(1,4)],
                                                             by = c("feature"="molecule"))) %>%
                mutate(cluster = as.numeric(as.character(cluster))) %>%
                mutate(cluster = ifelse(is.na(cluster), 0, cluster))

            #--5. compute silhouette
            sil <- silhouette(dmatrix, tmp.cluster$cluster)

            #--6. store
            result[result.index, "comp"] <- comp
            pos.res <-  sil$average.cluster  %>% 
                dplyr::filter(cluster == comp) %>% dplyr::pull(silhouette.coef)
            result[result.index, "pos"] <- ifelse(length(pos.res) == 0, NA, pos.res)
            neg.res <-  sil$average.cluster  %>%
                dplyr::filter(cluster == -comp) %>% dplyr::pull(silhouette.coef)
            result[result.index, "neg"] <- ifelse(length(neg.res) == 0, NA, neg.res)
            result.index <- result.index + 1
        }
    }
    result <- list("silhouette" = result)
    result[["ncomp"]] <- ncomp 
    result[["test.keepX"]] <- test.list.keepX
    result[["test.keepY"]] <- test.keepY
    result[["block"]] <- names(list.keepX.keepY)
    class(result) <- "block.spls.tune.silhouette"
    
    #-- 7. choice.keepX / choice.keepY
    result[["slopes"]] <- tune.silhouette.get_slopes(result)
    tmp <- tune.silhouette.get_choice_keepX(result) %>% 
        do.call(what = "rbind")
    if(!is.null(Y)){ # choice keepY
        result[["choice.keepY"]] <- tmp$Y
        tmp <- dplyr::select(tmp, -Y)
    }
    result[["choice.keepX"]] <- as.list(tmp)
    return(result)
}

#' @importFrom purrr imap set_names map_dfr
#' @importFrom dplyr filter mutate group_by summarise left_join n
#' @importFrom stringr str_split
tune.silhouette.get_slopes <- function(object){
    stopifnot(class(object) %in% c("block.spls.tune.silhouette", "spls.tune.silhouette", "spca.tune.silhouette"))
    # tune.silhouette is a data.frame (comp, X, Y, bock ..., pos, neg)
    # tune.silhouette <- tune.block.spls$silhouette
    
    block <- object$block
    ncomp <- object$ncomp
    
    if(is(object, "block.spls.tune.silhouette")){
        coord <- object$test.keepX
        if(!is.null(object[["test.keepY"]])){
            coord[["Y"]] <- object$test.keepY
            coord <- lapply(coord, sort)
        }
    }else if(is(object, "spls.tune.silhouette")){
        coord <- list(X= sort(object$test.keepX), 
                      Y= sort(object$test.keepY))
    }else if(is(object, "spca.tune.silhouette")){
        coord <- list(X= sort(object$test.keepX))
    }

    # get all points 
    all.points <- unique(object$silhouette[block])
    
    # define neighbours
    neighbourhood <- map_dfr(1:nrow(all.points), ~tune.silhouette.get_neighbours(coord = coord, all.points[.x,,drop=FALSE]))
    
    # extract pos and neg and set it as named list for better performance
    split_by_comp <- split(object$silhouette, f= as.factor(object$silhouette$comp))
    names_split_by_comp <- lapply(split_by_comp, 
                                  function(x) as.vector(apply(x[block], 1, function(y) paste(y, collapse = "_"))))
    
    POS <- purrr::imap(split_by_comp, ~ purrr::set_names(.x[["pos"]], names_split_by_comp[[.y]]))
    NEG <- purrr::imap(split_by_comp, ~ purrr::set_names(.x[["neg"]], names_split_by_comp[[.y]]))
    
    # compute slope (pos / neg by comp)
    slopes <- purrr::map_dfr(as.character(1:ncomp), ~{
        cbind(neighbourhood,
              "origin.pos" = POS[[.x]][neighbourhood$origin],
              "destination.pos" = POS[[.x]][neighbourhood$destination],
              "origin.neg" = NEG[[.x]][neighbourhood$origin],
              "destination.neg" = NEG[[.x]][neighbourhood$destination],
              "comp" = as.numeric(.x))})
    
    slopes <- slopes %>% dplyr::filter(origin!=destination) %>%
        dplyr::mutate("slope.pos" = tune.silhouette.get_slopes_coef(x1 = lapply(stringr::str_split(.$origin, "_"), as.numeric),
                                                             x2 = lapply(stringr::str_split(.$destination, "_"), as.numeric),
                                                             y1 = .$origin.pos,
                                                             y2 = .$destination.pos)) %>%
        dplyr::mutate("slope.neg" = tune.silhouette.get_slopes_coef(x1 = lapply(stringr::str_split(.$origin, "_"), as.numeric),
                                                             x2 = lapply(stringr::str_split(.$destination, "_"), as.numeric),
                                                             y1 = .$origin.neg,
                                                             y2 = .$destination.neg)) %>%
        dplyr::mutate("distance_from_origin" = tune.silhouette.distance_from_origin(x1 = lapply(stringr::str_split(.$origin, "_"), as.numeric)))
    
    # cumpute SD by comp and direction
    SD <- slopes %>% 
        dplyr::group_by(comp, direction) %>% 
        dplyr::summarise(sd.pos = sd.new(slope.pos, na.rm = TRUE),
                         sd.neg = sd.new(slope.neg, na.rm = TRUE),
                         mean.pos = mean(slope.pos, na.rm = TRUE),
                         mean.neg = mean(slope.neg, na.rm = TRUE), N=dplyr::n())
    
    # add Pval for signif slopes
    slopes <- slopes %>% dplyr::left_join(SD, by = c("direction", "comp")) %>%
        # pos
        dplyr::mutate(Z_score.pos = ifelse(.$sd.pos == 0,0,(.$slope.pos - .$mean.pos)/.$sd.pos)) %>%
        dplyr::mutate(Pval.pos = ifelse(Z_score.pos >= 0,1-pnorm(Z_score.pos), pnorm(Z_score.pos))) %>%
        # neg
        dplyr::mutate(Z_score.neg = ifelse(.$sd.neg == 0,0,(.$slope.neg - .$mean.neg)/.$sd.neg)) %>%
        dplyr::mutate(Pval.neg = ifelse(Z_score.neg >= 0,1-pnorm(Z_score.neg), pnorm(Z_score.neg)))
    return(slopes)
}

#' @importFrom purrr map2
tune.silhouette.get_neighbours <- function(coord, point){
    # return forward neighbourhood of a point
    
    # valid point in coord
    stopifnot(all(names(coord) == names(point)))
    
    neighbour.max <- point
    index <- purrr::map2(coord, point, ~which(.x == .y))

    # get max neighbours
    for(i in names(index)){
        if(!(length(coord[[i]]) == index[[i]])){ 
            neighbour.max[[i]] <- coord[[i]][index[[i]]+1]
        }
    } # can be simplified
    
    # get all close possible neighbours
    neighbourhood <- as.list(rbind(as.data.frame(point), as.data.frame(neighbour.max))) %>%
        lapply(unique) %>%
    expand.grid(stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
    
    # get direction for standard deviation computation
    direction <- vector(mode = "numeric", length = nrow(neighbourhood))
    for(i in 1:nrow(neighbourhood)){
        direction[i] <- purrr::map2(coord, neighbourhood[i,], ~which(.x == .y)) %>%
            purrr::map2(index, ~{.x-.y}) %>%
            paste(collapse = "_")
    }
    
    neighbourhood["direction"] <- direction
    neighbourhood["origin"] <- paste(point, collapse = "_")
    neighbourhood["destination"] <- apply(neighbourhood[names(coord)], 1, 
                                          function(x) paste(x, collapse = "_"))
    return(neighbourhood)
}

#' @importFrom purrr map_dbl
tune.silhouette.get_slopes_coef <- function(x1,x2,y1,y2){
    stopifnot(length(x1) == length(x2))
    stopifnot(length(x2) == length(y1))
    stopifnot(length(y1) == length(y2))
    
    res <- purrr::map_dbl(seq_along(x1), ~{
        euc.dist <- sqrt(sum((x1[[.x]] - x2[[.x]])^2));
        ifelse(euc.dist == 0, 0, (y2[.x] - y1[.x])/euc.dist)})
    return(res)
}

#' @importFrom purrr map_dbl
tune.silhouette.distance_from_origin <- function(x1){
    euc.dist <- purrr::map_dbl(seq_along(x1), ~{ sqrt(sum((x1[[.x]])^2))})
    return(euc.dist)
}

#' @importFrom dplyr select summarise left_join
#' @importFrom tidyr gather
#' @importFrom magrittr %>%
tune.silhouette.get_choice_keepX <- function(tune.block.spls){
    # from slopes, keep useful columns and remove NAs.
    slopes <- tune.block.spls$slopes %>% na.omit()
    tmp <- slopes %>%
        dplyr::select(c(tune.block.spls$block, comp, direction, Pval.pos, Pval.neg, distance_from_origin)) %>%
        tidyr::gather(Pval.dir, Pval.value, -c(tune.block.spls$block, comp, direction, distance_from_origin)) 
    
    # for each comp, arrange by Pvalue and distance from origin and get first result    
    MIN <- split(tmp, f=tmp$comp) %>% 
        lapply(function(x) x %>% dplyr::arrange(Pval.value, distance_from_origin) %>% .[1,] %>%
                   dplyr::select(tune.block.spls$block))
    
    return(MIN)
}
