# Jasmin Straube, Queensland Facility of Advanced Bioinformatics
#
# NOTE to myself: Kidney have been integrated to timeOmics after lmms removal from CRAN
# to remove WARNINGS, hidden variable ".Random.seed" have been removed and data have been re-saved.
# Raw data are still available: https://github.com/cran/lmms/blob/master/data/kidneySimTimeGroup.RData


#' Kidney Simulation Data
#' 
#' This data set contains the simulated expression of 140 proteins in 40 samples
#' from either group 1 or group 2 measured on the time points 0, 0.5, 1, 2, 3, 4. 
#' 
#' @usage 
#' data(kidneySimTimeGroup)
#' 
#' @description 
#' A list containing the following components:
#' @format
#' \describe{
#'    \item{\code{data}}{data matrix with 192 rows and 140 columns. Each row represents 
#'        an experimental sample, and each column a single protein.}
#'    \item{\code{time}}{a numeric vector containing the time points on which each sample is measured}
#'    \item{\code{sampleID}}{a character vector containing the biological replicate information of each sample }
#'    \item{\code{group}}{a character vector indicating the group of each sample}
#'    }
#'
#' @details
#'    This simulated data set consists of 40 samples and 140 proteins and was based on an
#'    the existing study from Freue \emph{et al.} (2010). Samples were measured on maximum 6 time points: 0, 0.5, 1, 2, 3, 4. Some samples have missing time points. 50 molecules were randomly selected to introduce a fold change of log(2). 
#'
#'
#' @source
#'    The Kidney Simulation Data is based  on the the paper of Freue \emph{et al.} (2010).
#'
#'
#' @references 
#'    Freue, G. V. C. et al. (2010). Proteomic signatures in plasma during early acute renal allograft rejection. \emph{Molecular & cellular proteomics}, \bold{9}, 1954-67.
#'
#'
#' @keywords 
#' datasets
"kidneySimTimeGroup"

