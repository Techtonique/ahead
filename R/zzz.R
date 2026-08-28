# zzz.R

#' Check that a Suggested package is installed, and error clearly if not
#'
#' @param pkg Name of the package required
#' @param fun Name of the calling function, for a clearer error message (optional)
#' @noRd
check_suggested <- function(pkg, fun = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      sprintf(
        "Package '%s' is required%s. Install it with install.packages('%s').",
        pkg,
        if (!is.null(fun)) paste0(" for ", fun) else "",
        pkg
      ),
      call. = FALSE
    )
  }
}
