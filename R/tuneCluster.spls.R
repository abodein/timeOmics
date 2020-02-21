#' tuneCluster.spls
#'
#' @examples
#' demo <- suppressMessages(get_demo_cluster())
#' X <- demo$X
#' Y <- demo$Y
#' tune.spls <- tuneCluster.spls(X, Y, ncomp = 2, test.keepX = c(5,10,15,20), test.keepY <- c(2,4,6))
#'
#' @export
#' @import mixOmics
#' @importFrom dplyr left_join
#' @importFrom dplyr mutate
#' @importFrom dplyr filter

tuneCluster.spls <- function(X, Y, ncomp = 2, test.keepX = rep(ncol(X), ncomp),
                             test.keepY = rep(ncol(Y), ncomp), ...){
    #-- checking input parameters ---------------------------------------------#
    #--------------------------------------------------------------------------#

    #-- X
    X <- validate.matrix.X(X)
    Y <- validate.matrix.Y(Y)

    #-- ncomp
    ncomp <- validate.ncomp(ncomp, list(X,Y))

    #-- test.keepX
    test.keepX <- validate.test.keepX(test.keepX = test.keepX, X = X)

    #-- test.keepY
    test.keepY <- validate.test.keepY(test.keepY = test.keepY, Y = Y)

    list.keepX.keepY <- list("keepX" = test.keepX, "keepY" = test.keepY) %>%
        expand.grid(stringsAsFactors = F, KEEP.OUT.ATTRS = F)


    #-- launch tuning  --------------------------------------------------------#
    #--------------------------------------------------------------------------#

    #-- 0. set output object
    result <- as.data.frame(matrix(ncol = 5, nrow = nrow(list.keepX.keepY)*ncomp))
    colnames(result) <- c("comp", "X", "Y", "pos", "neg")
    result.index <- 0

    #--1. compute dissimilarity matrix for silhouette coef. (once and for all)
    all_data <- cbind(X, Y)
    dmatrix <- dmatrix.spearman.dissimilarity(all_data)

    cluster <- as.data.frame(list("feature" = rownames(dmatrix),
                                  "block" = c(rep("X", ncol(X)), rep("Y", ncol(Y)))))

    #--2. tuning
    for(comp in 1:ncomp){
        for(index.list.kX.kY in 1:nrow(list.keepX.keepY)){
            # foreach comp, keepX and keepY of other comp is set to minimum
            kX <- rep(min(test.keepX), ncomp)
            kY <- rep(min(test.keepY), ncomp)
            kX[comp] <- list.keepX.keepY[index.list.kX.kY,"keepX"]
            kY[comp] <- list.keepX.keepY[index.list.kX.kY,"keepY"]

            #--3. run spls
            spls.res <- mixOmics::spls(X = X, Y = Y, ncomp =  ncomp, keepX = kX, keepY = kY, ...)

            #--4. extract clusters
            tmp.cluster <- getCluster(spls.res)
            tmp.cluster <- suppressWarnings(dplyr::left_join(cluster, tmp.cluster[c(1,4)],
                                                             by = c("feature"="molecule"))) %>%
                dplyr::mutate(cluster = as.numeric(as.character(cluster))) %>%
                dplyr::mutate(cluster = ifelse(is.na(cluster), 0, cluster))

            #--5. compute silhouette
            sil <- silhouette(dmatrix, tmp.cluster$cluster)

            #--6. store
            result.index <- result.index + 1
            result[result.index, "comp"] <- comp
            result[result.index, "X"] <- kX[comp]
            result[result.index, "Y"] <- kY[comp]
            result[result.index, "pos"] <- sil$average.cluster  %>%
                dplyr::filter(cluster == comp) %>% pull(silhouette.coef)
            result[result.index, "neg"] <- sil$average.cluster  %>%
                dplyr::filter(cluster == -comp) %>% pull(silhouette.coef)
        }
    }
    result <- list("silhouette" = result)
    result[["ncomp"]] <- ncomp
    result[["test.keepX"]] <- test.keepX
    result[["test.keepY"]] <- test.keepY
    result[["block"]] <- c("X", "Y")
    class(result) <- "spls.tune.silhouette"
    return(result)
}
