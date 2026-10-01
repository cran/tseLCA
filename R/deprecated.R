# tseLCA/R/deprecated.R
#
# Deprecation of the tseLCA 1.x interface. The 1.x functions keep working;
# each warns once per R session when called directly by the user (calls from
# within the package, which still uses them internally, do not warn). Set
# options(tseLCA.warn.deprecated = FALSE) to silence the warnings.

.tse_state <- new.env(parent = emptyenv())

#' Warn once per session that a 1.x function is deprecated
#' @noRd
.tse_deprecated <- function(what, replacement) {
  if (!isTRUE(getOption("tseLCA.warn.deprecated", TRUE))) {
    return(invisible(FALSE))
  }
  warned <- .tse_state$warned
  if (what %in% warned) {
    return(invisible(FALSE))
  }
  .tse_state$warned <- c(warned, what)
  warning(
    sprintf(
      "%s is deprecated as of tseLCA 2.0.0; use %s. %s",
      what,
      replacement,
      "(Shown once per session; see NEWS for the new interface.)"
    ),
    call. = FALSE
  )
  invisible(TRUE)
}

#' Deprecation warning for an exported 1.x helper, unless it was called from
#' within the package
#'
#' `frame` is the calling frame of the deprecated function (its
#' `parent.frame()`).
#' @noRd
.tse_deprecated_external <- function(what, replacement, frame) {
  if (identical(topenv(frame), topenv(environment(.tse_deprecated)))) {
    return(invisible(FALSE))
  }
  .tse_deprecated(what, replacement)
}
