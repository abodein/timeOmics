#' tuneCluster.block.spls
#'
#' @examples
#' demo <- suppressMessages(get_demo_cluster())
#' X <- list(X = demo$X, Z = demo$Z)
#' Y <- demo$Y
#' test.list.keepX <- list("X" = c(5,10,15,20), "Z" = c(2,4,6,8))
#' test.keepY <- c(2:5)
#' tune.block.spls <- tuneCluster.block.spls(X = X, Y = Y, test.list.keepX = test.list.keepX, test.keepY = test.keepY, mode = "canonical")

#' @import mixOmics
#' @importFrom dplyr left_join
#' @importFrom dplyr mutate
#' @importFrom dplyr filter
#' @export
tuneCluster.block.spls <- function(X, Y = NULL, ncomp = 2, test.list.keepX = rep(ncol(X), ncomp),
                                   test.keepY = NULL, indY = NULL, ...){
    #-- checking input parameters ---------------------------------------------#
    #--------------------------------------------------------------------------#

    #-- X
    X <- validate.list.matrix.X(X)

    #-- Y
    if(!is.null(Y)){
        Y <- validate.matrix.X(Y)
    } else {
        indY <- validate.indY(indY, X)
    }

    #-- ncomp
    ncomp <- validate.ncomp(ncomp, X)

    #-- keepX
    test.list.keepX <- validate.test.list.keepX(test.list.keepX, X, ncomp)

    #-- keepY
    test.keepY <- validate.test.keepX(test.keepY, Y, ncomp)

    list.keepX.keepY <- test.list.keepX
    list.keepX.keepY$Y <- test.keepY
    list.keepX.keepY <- expand.grid(list.keepX.keepY, stringsAsFactors = F,
                                    KEEP.OUT.ATTRS = F)


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
    result[["slopes"]] <- tune.silhouette.get_slopes(result)
    return(result)
}

#' @importFrom purrr imap set_names map_dfr
#' @importFrom dplyr filter mutate group_by summarise left_join
#' @importFrom stringr str_split
tune.silhouette.get_slopes <- function(object){
    stopifnot(class(object) %in% c("block.spls.tune.silhouette"))
    # tune.silhouette is a data.frame (comp, X, Y, bock ..., pos, neg)
    # tune.silhouette <- tune.block.spls$silhouette
    
    block <- object$block
    coord <- object$test.keepX
    if(!is.null(object[["test.keepY"]])){
        coord[["Y"]] <- object$test.keepY
    }
    coord <- lapply(coord, sort)
    ncomp <- object$ncomp
    
    # get all points 
    all.points <- unique(object$silhouette[block])
    
    # define neighbours
    neighbourhood <- map_dfr(1:nrow(all.points), ~tune.silhouette.get_neighbours(coord = coord, all.points[.x,]))
    
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
        dplyr::mutate("distance_from_origin" = tune.silhouette.distance_from_origin(x = lapply(stringr::str_split(.$origin, "_"), as.numeric)))
    
    # cumpute SD by comp and direction
    SD <- slopes %>% 
        dplyr::group_by(comp, direction) %>% 
        dplyr::summarise(sd.pos = sd(slope.pos),
                         sd.neg = sd(slope.neg),
                         mean.pos = mean(slope.pos),
                         mean.neg = mean(slope.neg))
    
    # add Pval for signif slopes
    slopes <- slopes %>% dplyr::left_join(SD, by = c("direction", "comp")) %>%
        # pos
        dplyr::mutate(Z_score.pos = (.$slope.pos - .$mean.pos)/.$sd.pos) %>%
        dplyr::mutate(Pval.pos = ifelse(Z_score.pos >= 0,1-pnorm(Z_score.pos), pnorm(Z_score.pos))) %>%
        # neg
        dplyr::mutate(Z_score.neg = (.$slope.neg - .$mean.neg)/.$sd.neg) %>%
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
    expand.grid(stringsAsFactors = F, KEEP.OUT.ATTRS = F)
    
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

tune.silhouette.choice_keepX <- function(slopes, alpha){
    # filter with pvalue on slope coef.
    slope.signif <- slopes %>% dplyr::filter(Pval.pos < alpha | Pval.neg <alpha) 
}

#' @importFrom dplyr mutate select filter
#' @importFrom tidyr gather
#' @import ggplot2
plot.block.spls.tune.silhouette <- function(object, pvalue = 0.05){
    data.gather <- object$silhouette %>%  
        dplyr::mutate(dim1 = pull(object$silhouette[block[1]])) %>% 
        dplyr::select(-c(object$block[1])) %>%
        tidyr::gather(dim2, value, -c(comp, pos, neg, dim1)) %>%
        dplyr::filter(comp == 1) %>% group_by(dim2,value)
    
    data.gather <- object$silhouette %>% mutate(dim2 = paste(Z,Y, sep = "_"))
    ggplot(data.gather, aes(x = X, y = pos, group = dim2)) + geom_line() + facet_wrap(~comp)
    
    data.gather <-  object$silhouette %>% gather(dim, value, -c(comp, pos, neg))
    
    ggplot(data.gather, aes(x = dim1, y = pos)) + geom_line(aes(col = dim2))
    
    plot.df <- map_dfr(object$block, ~{object$silhouette %>%
            unite("dim1", .x, remove = FALSE) %>%
            mutate(dim1 = as.numeric(dim1)) %>%
            tidyr::unite("dim2", one_of(object$block[which(object$block != .x)]), remove = FALSE) %>%
            mutate(dim2 = paste(.x, dim2, sep = "_")) %>%
            mutate(block = .x) %>%
            gather(silhouette, value, c(pos,neg))})
           
    ggplot(plot.df %>% filter( dim2 == "X_2_2") #Z == 2, Y == 2, comp==1, silhouette =="neg")
           ,aes(x = dim1, y = value, group = dim2, color = silhouette)) + geom_line() + facet_grid(comp~block, scales = "free_x")
    ggplot(plot.df, aes(x = dim1, y = value, group = dim2, color = silhouette)) + geom_line() + facet_grid(comp~block, scales = "free_x")
}

