#' Time-Course Simulated data
#' 
#' A list of dataset containing simulated data to be used in the vignette, generated as follow.
#' Twenty reference time profiles were generated on nine equally spaced time points and assigned to four clusters (five profiles each). 
#' These ground truth profiles were then used to simulate new profiles (5 each).
#' Finally, we modelled expression of new sampled profiles as a function of time.
#' 
#' @format a list of data.frame
#' \describe{
#'   \item{rawdata}{data.frame, reference profiles}
#'   \item{sim}{data.frame, sampled profiles from reference data}
#'   \item{modelled}{data.frame, modelled data}
#'  }
#' 
"timeOmics.simdata"