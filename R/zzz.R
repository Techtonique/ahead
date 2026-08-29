#' Check that a suggested package is installed
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
check_suggested <- function(pkg) {
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
