#' Check that a suggested package is installed, if not, install it 
#'
#' Internal helper used to guard functions that depend on packages listed
#' under \code{Suggests} rather than \code{Imports}. Throws an informative
#' error if the package is not installed, so users only need to install it
#' when they actually use a function that requires it. The install command
#' in the error message includes the Techtonique r-universe repository as a
#' fallback, so it works for both CRAN packages and Techtonique-hosted
#' packages that aren't on CRAN.
#'
#' @param pkg Character string naming the package to check (e.g. \code{"forecast"}).
#'
#' @return Invisibly returns \code{TRUE} if the package is available.
#'   Called for its side effect (throwing an error) otherwise.
#'
#' @noRd
check_suggested <- function(pkg, ask = interactive()) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    return(invisible(TRUE))
  }

  do_install <- TRUE
  if (ask) {
    do_install <- utils::askYesNo(
      sprintf("Package '%s' is required but not installed. Install it now?", pkg)
    )
    do_install <- isTRUE(do_install)
  }

  if (do_install) {
    utils::install.packages(
      pkg,
      repos = c("https://techtonique.r-universe.dev", "https://cloud.r-project.org")
    )
  }

  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      sprintf(
        "Package '%s' is required. Install it with install.packages('%s', repos = c('https://techtonique.r-universe.dev', 'https://cloud.r-project.org')).",
        pkg, pkg
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}