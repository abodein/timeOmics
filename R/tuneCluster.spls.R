#' Feature Selection Optimization for sPLS method 
#' 
#' This function identify the number of feautures to keep per component and thus by cluster in \code{mixOmics::spls} by optimizing the silhouette coefficient, which assesses the quality of clustering.
#' 
#' @param X numeric matrix (or data.frame) with features in columns and samples in rows
#' @param Y numeric matrix (or data.frame) with features in columns and samples in rows (same rows as \code{X})
#' @param ncomp integer, number of component to include in the model
#' @param test.keepX vector of integer containing the different value of keepX to test for block \code{X}.
#' @param test.keepY vector of integer containing the different value of keepY to test for block \code{Y}.
#' @param ... other parameters to be included in the spls model (see \code{mixOmics::spls})
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
#' \code{\link[mixOmics]{spls}}, \code{\link[timeOmics]{getCluster}}, \code{\link[timeOmics]{plotLong}}
#'
#' @examples
#' demo <- suppressWarnings(get_demo_cluster())
#' X <- demo$X
#' Y <- demo$Y
#' 
#' # tuning
#' tune.spls <- tuneCluster.spls(X, Y, ncomp= 2, test.keepX= c(5,10,15,20), test.keepY= c(2,4,6))
#' keepX <- tune.spls$choice.keepX
#' keepY <- tune.spls$choice.keepY
#' 
#' # final model
#' spls.res <- mixOmics::spls(X, Y, ncomp= 2, keepX= keepX, keepY= keepY)
#' 
#' # get clusters and plot longitudinal profile by cluster
#' spls.cluster <- getCluster(spls.res)
#' plotLong(spls.res)
#' 
#'
#' @export
#' @import mixOmics
#' @importFrom dplyr left_join
#' @importFrom dplyr mutate
#' @importFrom dplyr filter

tuneCluster.spls <- function(X, Y, ncomp = 2, test.keepX = rep(ncol(X), ncomp),
                             test.keepY = rep(ncol(Y), ncomp), ...){
    #-- checking input parameters ------------  ---------------------------------#
    #--------------------------------------------------------------------------#

    #-- X
    X <- validate_matrix_X(X)
    Y <- validate_matrix_Y(Y)

    #-- ncomp
    ncomp <- validate_ncomp(ncomp, list(X,Y))

    #-- test.keepX
    test.keepX <- validate_test_keepX(test.keepX = test.keepX, X = X)

    #-- test.keepY
    test.keepY <- validate_test_keepY(test.keepY = test.keepY, Y = Y)

    list.keepX.keepY <- list("keepX" = test.keepX, "keepY" = test.keepY) %>%
        expand.grid(stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)


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
            # result[result.index, "pos"] <- sil$average.cluster  %>%
            #     dplyr::filter(cluster == comp) %>% pull(silhouette.coef)
            # result[result.index, "neg"] <- sil$average.cluster  %>%
            #     dplyr::filter(cluster == -comp) %>% pull(silhouette.coef)
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
    result[["test.keepY"]] <- test.keepY
    result[["block"]] <- c("X", "Y")
    class(result) <- "spls.tune.silhouette"
    
    #-- 7. choice.keepX / choice.keepY
    result[["slopes"]] <- tune.silhouette.get_slopes(result)
    tmp <- tune.silhouette.get_choice_keepX(result) # choice keepX/keepY
    result[["choice.keepX"]] <- unlist(lapply(tmp, function(x) x$X))
    result[["choice.keepY"]] <- unlist(lapply(tmp, function(x) x$Y))
    return(result)
}
