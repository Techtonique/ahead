.onLoad <- function(libname, pkgname){
 if (!("ForecastComb" %in% rownames(utils::installed.packages()))) 
 {
    utils::install.packages("ForecastComb", repos="https://cloud.r-project.org")
 }
}