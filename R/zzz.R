.onAttach <- function(libname, pkgname) {
    invisible(suppressPackageStartupMessages(
        sapply(c("tibble", "purrr", "dplyr", "tidyr", "ggplot2", "mixOmics"),
               requireNamespace, quietly = TRUE)
    ))
}