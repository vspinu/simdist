##' @include distances.R
NULL

##' Similarities for sparse and dense 2d representations
##'
##' Similarity counterparts of \code{\link{distances}}, which see. Generally
##' similarities are computed from their distance counterparts through
##' \code{s=1-d} or \code{s=1/(1+d)}, but this need not always be the case. See
##' the source code for the exact computation.
##' 
##' @inheritParams dist_cosine
##' @rdname similarities
##' @export
sim_cosine <- function(x, y = NULL, ptrans = "l2", strans = NULL,
                       by = c("primary", "secondary", "row", "column"),
                       ref_names = NULL) {
    1 - generic_dist(x, y, ptrans, strans, by, ref_names, dist_type = "COSINE")
}

##' @rdname similarities
##' @export
sim_euclidean <- function(x, y = NULL, ptrans = NULL, strans = NULL,
                          by = c("primary", "secondary", "row", "column"),
                          ref_names = NULL) {
    1/(1 + generic_dist(x, y, ptrans, strans, by, ref_names, dist_type = "EUCLIDEAN"))
}


##' Aggregation similarities for sparse 3d representations 
##'
##' These are similarity counterparts of \code{aggregated_distances}, which
##' see. Generally similarities are computed from their distance counterparts
##' through \code{s=1-d} or \code{s=1/(1+d)}, but this need not always be the
##' case. See the source code for the exact computation.
##'
##' @inheritParams adist_centroid
##' @rdname aggregated_similarities
##' @export
asim_centroid <- function(x, y = NULL, vecs, ptrans = NULL, strans = NULL,
                          by = c("primary", "secondary", "row", "column"),
                          precompute = TRUE, dist_type = "cosine") {
    dist <- generic_aggr_dist(x, y, vecs, ptrans, strans, by = by, precompute = precompute,
                              aggr_type = "CENTROID", dist_type = dist_type)
    if (dist_type == "cosine")
        1 - dist
    else
        1/(1 + dist)
}


##' @rdname aggregated_similarities
##' @export
asim_semantic_min_sum <- function(x, y = NULL, vecs, ptrans = NULL, strans = NULL,
                                  by = c("primary", "secondary", "row", "column"),
                                  precompute = TRUE, dist_type = "cosine") {
    dist <- generic_aggr_dist(x, y, vecs, ptrans, strans, by = by, precompute = precompute,
                              aggr_type = "SEM_MIN_SUM", dist_type = dist_type)
    if (dist_type == "cosine")
        2 - dist
    else
        1/(1 + dist)
}


##' @rdname aggregated_similarities
##' @export
asim_semantic_min_max <- function(x, y = NULL, vecs, ptrans = NULL, strans = NULL,
                                  by = c("primary", "secondary", "row", "column"),
                                  precompute = TRUE, dist_type = "cosine") {
    dist <- generic_aggr_dist(x, y, vecs, ptrans, strans, by = by, precompute = precompute,
                              aggr_type = "SEM_MIN_MAX", dist_type = dist_type)
    if (dist_type == "cosine")
        2 - dist
    else
        1/(1 + dist)
}

##' @rdname aggregated_similarities
##' @export
asim_rwmd <- asim_semantic_min_max
