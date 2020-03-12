.onAttach <- function(libname, pkgname) {
    invisible(suppressPackageStartupMessages(
        for( i in c("tibble", "purrr", "dplyr", "tidyr", "ggplot2", "mixOmics")){
            requireNamespace(i, quietly = TRUE)
        }
    ))
}