##' High performance computation of similarities and distances for sparse and
##' dense objects.
##'
##' The package contains efficient parallel functions for computation of
##' similarity and distance metrics on various sparse and dense
##' representations. Canonical applications of these functions are natural
##' language processing and recommender systems.
##'
##' \code{Simdist} package uses a higher level abstraction for 2d sparse
##' representation than the standard sparse matrices software. For every
##' supported 2d representation \code{primary} and \code{secondary} dimension of
##' variation of the measurement are defined. Every function in this package
##' acts either on primary or secondary dimension. The primary reason for
##' primary/secondary division is computational - computing along primary
##' dimension is usually more efficient than along the secondary dimension. Even
##' for dense matrices the "mental model" used in the package is that of nested
##' lists - higher order grouping (i.e. document) is stored as entries along
##' primary dimension and inner elements (i.e. terms) are stored across
##' secondary dimension.
##'
##' @section Supported representations:
##' 
##' The supported 2d representations and primary-secondary dimensions are as
##' follows:
##'
##' \describe{
##'
##' \item{matrix}{primary - columns, secondary - rows}
##'
##' \item{\code{Matrix::dgCMatrix}}{primary - columns, secondary - rows}
##'
##' \item{\code{Matrix::dgRMatrix}}{primary - rows, secondary - columns}
##'
##' \item{\code{Matrix::dgTMatrix}}{primary - rows, secondary - columns}
##'
##' \item{\code{slam::simple_triplet_matrix}}{primary - rows, secondary -
##'   columns (not yet supported)}
##'
##' \item{data frames in primary-secondary-value (psv) format}{primary - first
##'   id column, secondary - second id column. Id and value columns could be
##'   explicitly marked with \code{\link{psv}} function.}
##'
##' \item{list of named numeric or character vectors}{primary - first list
##'   level, secondary - inner vector level (not yet implemented)}
##'
##' }
##'
##' To minimize the risk of logical errors due to mismatched dimensions only
##' distances across same-type objects are currently implemented.
##'
##' @section Non-conventional approach:
##' \itemize{
##' 
##' \item{Primary/secondary dimension distinction which allows treating all
##'   representations as two-level nested lists.}
##' 
##' \item{For named matrices,  secondary dimensions are matched by names, not
##'   positionaly. This means that even for matrices the size of the secondary
##'   dimension need not match. All rows in X not in Y will be considered
##'   missing (aka 0s) as if it were a sparse matrix.}
##'
##' \item{No normalization by default. All sim and dist functions accept
##'   normalization or scaling functions (transformers) which allow arbitrary
##'   transformation of the input matrices.}
##'
##' \item{Cosine similarity of a vector X with 0 vector is 0,  in contrast to
##'   \code{proxy} package where it's 1. This preserves coordinate-wise
##'   continuity in 0 and allows for a more efficient implementation.}
##'
##' }
##' 
##' @docType package
##' @name simdist-package
##' @aliases simdist
##' @import Matrix
##' @importFrom methods setClass setGeneric setMethod new show callGeneric
##' @importFrom Rcpp sourceCpp
##' @importFrom RcppParallel RcppParallelLibs
##' @importFrom fastmatch fmatch
##' @useDynLib simdist
##' @include transform.R
##' @include utils.R
##' @include S4utils.R
NULL


### Generic wrappers on objects of same input class

dist_same_mat_class <- function(x, y, ptrans, strans,
                                by = c("primary", "secondary", "row", "column"),
                                ref_names = NULL, pairwise, dist_fun) {
    ## TODO: if rownames are not given, enforce same matrix dimmension
    by <- match.arg(by)
    x <- to_primary_maybe(.norm(x, ptrans, strans), by)
    y <- to_primary_maybe(.norm(y, ptrans, strans), by)
    xsec_null <- is.null(secondary_names(x))
    ysec_null <- is.null(secondary_names(y))
    if (xsec_null && !ysec_null)  warning ("X has no secondary names, but Y does; discarding.")
    if (ysec_null && !xsec_null)  warning ("Y has no secondary names, but X does; discarding.")
    if (xsec_null || ysec_null) {
        if (secondary_size(x) != secondary_size(y)) {
            stop ("When matrices have no secondary names, secondary dimmension must match.")
        }
        out <- dist_fun(x, secondary_ix(x), y, secondary_ix(y), 0L, 0L)
    } else {
        xresort <- 0L
        yresort <- 0L
        if (is.null(ref_names)) {
            ## no nas here 
            ref_names <- c(secondary_names(x), secondary_names(y))
            xix <- secondary_ix(x)
            yresort <- 2L
        } else {
            xnnas <- !is.na(fastmatch::fmatch(secondary_names(x), ref_names))
            ynnas <- !is.na(fastmatch::fmatch(secondary_names(y), ref_names))
            if (!any(xnnas))
                stop("No names in X matched ref_names. Are you computing on the right dimension?")
            if (!any(ynnas))
                stop("No names in Y matched ref_names. Are you computing on the right dimension?")
            x <- secondary_subset(x, xnnas)
            y <- secondary_subset(y, ynnas)
            xix <- fastmatch::fmatch(secondary_names(x), ref_names)[secondary_ix(x) + 1L] - 1L
            ## subset copied x and y, re-sorting in place
            xresort <- 1L
            yresort <- 1L
        }
        yix <- fastmatch::fmatch(secondary_names(y), ref_names)[secondary_ix(y) + 1L] - 1L
        out <- dist_fun(x, xix, y, yix, xresort, yresort)
    }
    xnames <- primary_names(x)
    ynames <- primary_names(y)
    if (pairwise) {
        data.frame(x = if(is.null(xnames)) 1:primary_size(x) else xnames,
                   y = if(is.null(ynames)) 1:primary_size(y) else ynames,
                   val = c(out_mat))
    } else {
        if(!is.null(xnames) || !is.null(ynames))
            dimnames(out) <- list(xnames, ynames)
        out
    }
}

dist_df <- function(x, y, ptrans = NULL, strans = NULL,
                    by = c("primary", "secondary", "row", "column"),
                    ref_names = NULL, pairwise, dist_fun) {
    by <- match.arg(by)
    x <- to_primary_maybe(.norm(x, ptrans, strans), by)
    y <- to_primary_maybe(.norm(y, ptrans, strans), by)
    xsec <- secondary_ids(x)
    ysec <- secondary_ids(y)
    if (is.null(ref_names)) {
        ref_names <- sort(unique(c(secondary_names(x), secondary_names(y))))
        check_nas <- FALSE
    } else {
        check_nas <- TRUE
    }

    xix <-
        if (is.factor(xsec)) fastmatch::fmatch(levels(xsec), ref_names)[xsec]
        else fastmatch::fmatch(xsec, ref_names)
    yix <-
        if (is.factor(ysec)) fastmatch::fmatch(levels(ysec), ref_names)[ysec]
        else fastmatch::fmatch(ysec, ref_names)

    xnames <- primary_names(x)
    ynames <- primary_names(y)

    out_mat <-
        if (check_nas) {
            xnnas <- !is.na(xix)
            ynnas <- !is.na(yix)
            dist_fun(primary_ix(x)[xnnas], xix[xnnas] - 1L, sparse_vals(x)[xnnas],
                     primary_ix(y)[ynnas], yix[ynnas] - 1L, sparse_vals(y)[ynnas],
                     length(xnames), length(ynames))
        } else {
            dist_fun(primary_ix(x), xix - 1L, sparse_vals(x),
                     primary_ix(y), yix - 1L, sparse_vals(y),
                     length(xnames), length(ynames))
        }
    
    ## out <- data.frame(x = xnames[rep(1:length(xnames), each = length(ynames))],
    ##                   y = ynames[rep(1:length(ynames), times = length(ynames))],
    ##                   dist = c(out_mat))
    ## class(out) <- class(x)

    if (pairwise) {
        structure(list(x = xnames,
                       y = ynames,
                       val = c(out_mat)),
                  class = class(x))
    } else {
        dimnames(out_mat) <- list(xnames, ynames)
        out_mat
    }
}


 ### GENERIC DISTANCE

setGeneric("generic_dist", signature = c("x", "y"),
           function(x, y, ptrans = NULL, strans = NULL,
                    by = c("primary", "secondary", "row", "column"),
                    ref_names = NULL, pairwise = FALSE, ..., dist_type) {
               standardGeneric("generic_dist")
           })

.generic_2d_mat_dist <- function(x, y, ptrans, strans,
                                 by = c("primary", "secondary", "row", "column"),
                                 ref_names = NULL, pairwise = FALSE, ..., dist_type) {
    fun <- function(x, xix, y, yix, xreorder, yreorder) {
        self <- FALSE
        if (is.matrix(x)) {
            ## fixme: move this to C
            if (xreorder) {
                xord <- order(xix)
                xix <- xix[xord]
                x <- x[xord, ]
            }
            if (yreorder) {
                yord <- order(yix)
                yix <- yix[yord]
                y <- y[yord, ]
            }
            C_dense_dist(dist_type, xix, x, yix, y, pairwise, self)
        } else if (inherits(x, "dgTMatrix")) {
            C_triplet_dist(dist_type, x@i, xix, x@x, y@i, yix, y@x,
                           primary_size(x), primary_size(y),
                           pairwise, self,
                           xreorder, yreorder)
        } else if (inherits(x, "dsparseMatrix")) {
            ## x and y can be either RsparseMatrix or CsparseMatrix
            C_sparse_dist(dist_type,
                          xix, x@p, x@x, yix, y@p, y@x,
                          primary_size(x), primary_size(y),
                          pairwise, self,
                          xreorder, yreorder)
        } else stop("unsuported matrix type")
    }
    dist_same_mat_class(x, y, ptrans, strans, by, ref_names = ref_names,
                        pairwise = pairwise, dist_fun = fun)
}

setMethod("generic_dist", signature("ANY", "NULL"),
          function(x, y, ptrans = NULL, strans = NULL,
                   by = c("primary", "secondary", "row", "column"),
                   ref_names = NULL, pairwise = FALSE, ..., dist_type) {
              generic_dist(x, x, ptrans = ptrans, strans = strans, by = by,
                           ref_names = ref_names, pairwise = pairwise, ..., dist_type = dist_type)
          })
setMethod("generic_dist", signature("matrix", "matrix"), .generic_2d_mat_dist)
setMethod("generic_dist", signature("dgRMatrix", "dgRMatrix"), .generic_2d_mat_dist)
setMethod("generic_dist", signature("dgCMatrix", "dgCMatrix"), .generic_2d_mat_dist)
setMethod("generic_dist", signature("dgTMatrix", "dgTMatrix"), .generic_2d_mat_dist)
setMethod("generic_dist", signature("data.frame", "data.frame"), {
    function(x, y, ptrans = NULL, strans = NULL,
             by = c("primary", "secondary", "row", "column"),
             ref_names = NULL, pairwise = pairwise, ..., dist_type) {
        self <- FALSE
        fun <- function(xpix, xsix, xval, ypix, ysix, yval, nrows, ncols) {
            C_triplet_dist(dist_type,
                           xpix, xsix, xval,
                           ypix, ysix, yval,
                           nrows, ncols,
                           pairwise, self, 
                           1L, 1L)
        }
        dist_df(x, y, ptrans, strans, by, ref_names = ref_names, pairwise = pairwise, dist_fun = fun)
    }
})


 ### GENERIC AGGREGATED DISTANCE

setGeneric("generic_aggr_dist", signature = c("x", "y"),
           function(x, y = NULL, vecs, ptrans = NULL, strans = NULL,
                    by = c("primary", "secondary", "row", "column"),
                    precompute = TRUE, pairwise = FALSE, ...,
                    aggr_type, dist_type) {
               standardGeneric("generic_aggr_dist")
           })

.generic_aggr_mat_dist <- function(x, y = NULL, vecs, ptrans, strans,
                                   by = c("primary", "secondary", "row", "column"),
                                   precompute = TRUE, pairwise = FALSE, ...,
                                   aggr_type, dist_type) {
    self <- FALSE
    aggr_type <- toupper(aggr_type)
    dist_type <- toupper(dist_type)
    fun <- function(x, xix, y, yix, ...) {
        if (inherits(x, "dgTMatrix")) {
            C_triplet_aggr_dist(aggr_type, dist_type, vecs,
                                x@i, xix, x@x, y@i, yix, y@x,
                                primary_size(x), primary_size(y),
                                pairwise, self, precompute)
        } else if (inherits(x, "dsparseMatrix")) {
            C_sparse_aggr_dist(aggr_type, dist_type, vecs,
                               xix, x@p, x@x, yix, y@p, y@x,
                               primary_size(x), primary_size(y),
                               pairwise, self, precompute)
        } else stop("unsuported matrix type")
    }
    dist_same_mat_class(x, y, ptrans, strans, by, ref_names = primary_names(vecs),
                        pairwise = pairwise, dist_fun = fun)
}

setMethod("generic_aggr_dist", signature("ANY", "NULL"),
          function(x, y = NULL, vecs, ptrans = NULL, strans = NULL,
                   by = c("primary", "secondary", "row", "column"),
                   ref_names = NULL, ..., aggr_type, dist_type) {
              generic_aggr_dist(x, x, vecs, ptrans = ptrans, strans = strans, by = by,
                                ref_names = ref_names, ...,
                                aggr_type = aggr_type, dist_type = dist_type)
          })
setMethod("generic_aggr_dist", signature("dgRMatrix", "dgRMatrix"), .generic_aggr_mat_dist)
setMethod("generic_aggr_dist", signature("dgCMatrix", "dgCMatrix"), .generic_aggr_mat_dist)
setMethod("generic_aggr_dist", signature("dgTMatrix", "dgTMatrix"), .generic_aggr_mat_dist)

setMethod("generic_aggr_dist", signature("data.frame", "data.frame"), {
    function(x, y, vecs, ptrans = NULL, strans = NULL,
             by = c("primary", "secondary", "row", "column"),
             precompute = TRUE, pairwise = FALSE, ...,
             aggr_type, dist_type) {
        self <- FALSE
        dist_type <- toupper(dist_type)
        fun <- function(xpix, xsix, xval, ypix, ysix, yval, nrows, ncols) {
            C_triplet_aggr_dist(aggr_type, dist_type, vecs,
                                xpix, xsix, xval,
                                ypix, ysix, yval,
                                nrows, ncols,
                                pairwise, self, precompute)
        }
        dist_df(x, y, ptrans, strans, by, ref_names = primary_names(vecs),
                pairwise = pairwise, dist_fun = fun)
    }
})
