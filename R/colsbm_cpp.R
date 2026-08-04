##' Low-level C++ bindings for colSBM VEM
##'
##' These are minimal wrappers around the C++ scaffold. The real
##' implementations will provide full functionality later.
##'
##' @param A List of adjacency matrices (numeric matrices)
##' @param Q Integer number of blocks
##' @param max_step Integer max iterations for optimize
##' @param tol Numeric tolerance
##' @return An external pointer to the internal C++ object (for create),
##' and other functions return diagnostics.
##' @name colsbm_cpp
NULL

#' @export
colsbm_create <- function(A, Q) {
    .Call(`_colSBM_colsbm_create`, A, as.integer(Q))
}

#' @export
colsbm_optimize <- function(ptr, max_step = 100L, tol = 1e-6) {
    .Call(`_colSBM_colsbm_optimize`, ptr, as.integer(max_step), as.double(tol))
}

#' @export
colsbm_vbound <- function(ptr) {
    .Call(`_colSBM_colsbm_vbound`, ptr)
}

#' @export
colsbm_info <- function(ptr) {
    .Call(`_colSBM_colsbm_info`, ptr)
}
