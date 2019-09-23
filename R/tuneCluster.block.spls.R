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
                                   test.keepY = NULL, ...){
    #-- checking input parameters ---------------------------------------------#
    #--------------------------------------------------------------------------#

    #-- X
    X <- validate.list.matrix.X(X)

    #-- Y
    if(!is.null(Y)){
        Y <- validate.matrix.X(Y)
    } else {
        indY <- validate.indY(indY)
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
            result[result.index, "pos"] <- sil$average.cluster  %>%
                dplyr::filter(cluster == comp) %>% dplyr::pull(silhouette.coef)
            result[result.index, "neg"] <- sil$average.cluster  %>%
                dplyr::filter(cluster == -comp) %>% dplyr::pull(silhouette.coef)
            result.index <- result.index + 1
        }
    }
    result <- list("silhouette" = result)
    result[["ncomp"]] <- ncomp
    result[["test.keepX"]] <- test.list.keepX
    result[["test.keepY"]] <- test.keepY
    result[["block"]] <- names(list.keepX.keepY)
    class(result) <- "block.spls.tune.silhouette"
    return(result)
}

tune.silhouette.get_slopes <- function(tune.silhouette){
    # tune.silhouette is a data.frame (comp, X, Y, bock ..., pos, neg)
    # tune.silhouette <- tune.block.spls$silhouette
    block <- names(tune.silhouette)[!names(tune.silhouette) %in% c("comp", "pos", "neg")]
    tune.silhouette[block]
}



