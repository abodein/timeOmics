#' lmms.filter.lines
#' 
#' This function filters linear models with highly heterogeneous variability within residues.
#' From an "lmms" output, 2 parameters are tested:
#' 
#' * homo-sedasticity of the residues with a Breusch-Pagan test 
#' * low dispersion with a cutoff on the MSE (mean squared error)
#' 
#' @param data a data.frame used in the \code{lmms::lmmSpline} command
#' @param lmms.obj a \code{lmmspline} object
#' @param time a numeric vector containing the sample time point information.
#' @param homoskedasticity a logical whether or not to test for homoscedasticity with the Breusch-Pagan test.
#' @param MSE.filter whether or not to test for low dispersion with a cutoff on the MSE.
#' @param homoskedasticity.cutoff a numeric scalar between 0 and 1, p-value threshold for B-P test.
#' 
#' @return 
#' a list containing the following items
#' \item{filtering.summary}{a data.frame with the different tests per features (passed = TRUE, failed = FALSE)}
#' \item{to.keep}{features which passed all the tests}
#' \item{filtered}{the filtered data.frame}
#' 
#' @seealso 
#' \code{\link[lmms]{lmmSpline}}, \code{\link[lmtest]{bptest}}
#' 
#' @examples 
#' # data and lmms output
#' data(timeOmics.simdata)
#' data <- timeOmics.simdata$sim
#' lmms.output <- timeOmics.simdata$lmms.output
#' time <- timeOmics.simdata$time
#' 
#' # filter
#' filter.res <- lmms.filter.lines(data = data, lmms.obj = lmms.output, time = time)
#'
#' 
#' @importFrom lmtest bptest
#' @importFrom purrr is_scalar_vector set_names
#' @importFrom  dplyr select mutate right_join filter left_join
#' 
#' @export
lmms.filter.lines <- function(data, 
                              lmms.obj, 
                              time,
                              homoskedasticity = TRUE, 
                              MSE.filter = TRUE, 
                              homoskedasticity.cutoff = 0.05)
{
    #-- Check parameters
    #------------------------------------
    data <- validate.matrix.X(data)
    
    #-- lmms.obj
    if(!is(lmms.obj, "lmmspline")) stop("'lmms.obj' should be a 'lmms' object.")
    if(!(any(slot(lmms.obj, "modelsUsed") == 0))){
        res <- list()
        res[["to.keep"]] <- colnames(data)
        res[["filtered"]] <- data
        return(res)
    }
    if(purrr::is_empty(slot(lmms.obj, "models"))){
        stop("No models found in 'lmms.obj', please use 'keepModels = TRUE' in 'lmmSpline()'")
    }
    
    #-- time
    if(!is.numeric(time) || length(time) != nrow(data)){
        stop("'time' should be a numeric vector with the same length as 'nrow(data)'")
    }
    if(any(!(time %in% as.numeric(colnames(slot(lmms.obj, "predSpline")))))){
        stop("wrong time between 'lmms.obj', and 'time'")
    }
    
    #-- homoskedasticity  TRUE/FALSE
    homoskedasticity = return.true.false(x = homoskedasticity, default = TRUE)
    
    #-- MSE.filter
    MSE.filter = return.true.false(x = MSE.filter, default = TRUE)
    
    #-- homoskedasticity.cutoff
    if(!purrr::is_scalar_vector(homoskedasticity.cutoff) || !is.finite(homoskedasticity.cutoff) || !is.numeric(homoskedasticity.cutoff)){
        stop("'homoskedasticity.cutoff' should be a numeric between 0 and 1")
    }
    if(homoskedasticity.cutoff > 1 || homoskedasticity.cutoff < 0){
        stop("'homoskedasticity.cutoff' should be a numeric between 0 and 1")  
    }
    
    #-- Run filter
    #--------------------------------------------
    predSpline <- slot(lmms.obj,'predSpline') %>% t %>% as.data.frame()
    modelsUsed <- slot(lmms.obj,'modelsUsed')
    if(all(modelsUsed == 0)){
        MSE.filter <- FALSE
    }
    
    result <- colnames(predSpline) %>% as.data.frame() %>%
        purrr::set_names(c("feature")) %>%
        dplyr::mutate(modelsUsed = modelsUsed, 
                      feature = as.character(feature))
    
    # homoskedasticity : Breusch-Pagan test
    # ---------------------
    # WARNING : only for linear regression model (modelsUsed == 0)
    # * if pvalue < signif cutoff : heteroskedasticity
    # * if pvalue > signif cutoff : homooskedasticity  -> to keep
    if(homoskedasticity){
        models0 <- slot(lmms.obj, "models")[modelsUsed == 0]
        
        BP.res <- lapply(models0, lmtest::bptest) %>% 
            lapply(function(x) as.numeric(x$p.value)) %>% unlist() %>%
            as.data.frame() %>% set_names("BP.test") %>%
            dplyr::mutate(feature = colnames(predSpline)[modelsUsed == 0])
        
        result <- result %>% dplyr::left_join(BP.res, by = c("feature", "feature")) %>% 
            dplyr::mutate(BP.test = ifelse(is.na(BP.test), 1, BP.test)) %>% 
            # replace (model != 0) by  a p value of 1, must be homoskedastic if not, no spline
            dplyr::mutate(BP.test = (BP.test >=  homoskedasticity.cutoff)) # T/F
    }
    
    # applied filter based on max MSE for model != 0
    if (MSE.filter) {
        MSE.res <- get_MSE(data = data, predSpline = predSpline, time = time, modelsUsed = modelsUsed)
        MSE.cutoff <- MSE.res %>% dplyr::filter(modelsUsed != 0) %>% pull(MSE) %>% max
        result <- MSE.res %>% 
            dplyr::mutate(MSE.filter = (MSE <= MSE.cutoff)) %>% # TRUE
            dplyr::select(feature, MSE.filter) %>% 
            dplyr::right_join(result, by=c("feature"))
    }
    
    # summary
    OK <- sum(c(homoskedasticity, MSE.filter)) # nb test to pass
    
    to.keep <- tidyr::pivot_longer(data = result, names_to='test', 
                                   values_to='res', -c(feature, modelsUsed)) %>%
        dplyr::group_by(feature) %>% 
        summarise(val = sum(res)) %>%
        filter(val == OK) %>% # pass all tests
        pull(feature)
    
    
    res <- list()
    res[["filtering.summary"]] <- result
    res[["to.keep"]] <- to.keep
    res[["filtered"]] <- dplyr::select(predSpline, to.keep)
    
    return(res)
}

#' @importFrom dplyr mutate group_by left_join summarise
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
get_MSE <- function(data, predSpline, time, modelsUsed){
    
    X1 <- data %>% as.data.frame() %>% dplyr::mutate(time = time) %>%
        tidyr::pivot_longer(names_to = "feature", values_to = "Y_i", -time)
    
    X2 <- predSpline %>% t %>% as.data.frame() %>%
        tibble::rownames_to_column("feature") %>%
        dplyr::mutate(modelsUsed = modelsUsed) %>%
        tidyr::pivot_longer(names_to = "time", values_to = "Y_hat", -c(feature, modelsUsed)) %>%
        dplyr::mutate(time = as.numeric(time)) 
    
    MSE <- dplyr::left_join(X1,X2, by=c("time"="time", "feature"="feature")) %>%
        na.omit() %>%
        dplyr::mutate(error = (Y_i-Y_hat)^2) %>%
        dplyr::group_by(feature, modelsUsed) %>%
        dplyr::summarise(MSE = mean(error))
    
    return(MSE)
}
