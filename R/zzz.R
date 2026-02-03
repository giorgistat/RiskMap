# Package environment for caching
.dsgm_cache <- new.env(parent = emptyenv())

##' @keywords internal
.onLoad <- function(libname, pkgname) {
  invisible(NULL)
}
