.onLoad <- function(libname, pkgname){
    if (!requireNamespace("devtools", quietly = TRUE)) {
      utils::install.packages("devtools")
    }
    if (!requireNamespace("misc", quietly = TRUE)) {
    devtools::install_github("thierrymoudiki/misc")
    }
}