##' @include core.R
NULL


##' Distances for sparse 2d representations
##'
##' @param x,y sparse or dense objects supported by \code{simdist}. See
##'     \code{\link{simdist-package}} for description of supported types and
##'     related terminology.
##' @param by Dimension along which to perform distance computation. For all
##'     supported data structures computation along primary dimension is more or
##'     as efficient than along the secondary dimension.
##' @param ptrans,strans Primary and secondary transformations. Can be either a
##'     function, string or a numeric vector. When a function, it must take 3
##'     arguments - an object supported by \code{simdist} distance measures,
##'     names of a dimension ("primary" or "secondary") and logical na.rm. When
##'     a string it must name a transformation function without the semantic
##'     prefix (norm_, scale_ or trans_); see \code{link{transformations}}. When
##'     numeric vector it specifies weights to scale along the corresponding
##'     dimension.
##' @param ref_names reference names for secondary entries. By default both
##'     \code{x} and \code{y} secondary names are used. Secondary entries not in
##'     \code{ref_names} are discarded from \code{x} and \code{y}.
##' @rdname distances
##' @aliases distances
##' @export
dist_cosine <- function(x, y = NULL, ptrans = "l2", strans = NULL,
                        by = c("primary", "secondary", "row", "column"),
                        ref_names = NULL, pairwise = FALSE) {
    generic_dist(x, y, ptrans, strans, by, ref_names, pairwise, dist_type = "COSINE")
}

##' @rdname distances
##' @export
dist_euclidean <- function(x, y = NULL, ptrans = NULL, strans = NULL,
                           by = c("primary", "secondary", "row", "column"),
                           ref_names = NULL, pairwise = FALSE) {
    generic_dist(x, y, ptrans, strans, by, ref_names, pairwise, dist_type = "EUCLIDEAN")
}


##' Aggregation distances for sparse 3d representations
##'
##' Distances for situation when every entry on secondary dimension is
##' characterized by a numeric vector (embedding). In the example of
##' term-document matrix where document is a primary dimension, each term has a
##' numeric representation in a N-dimensional space. For user-movie rating,
##' vectors for movies can represent various movie characteristics. The
##' aggregation distances (adist for short) perform various aggregation steps of
##' these vectors
##'
##' \describe{
##'
##' \item{\code{centroid}}{Within each primary entry (document, user etc.) the
##'   vectors of secondary entries (terms, movies etc) are averaged element-wise
##'   and \code{dist_type} is applied on the resulting vectors.}
##'
##' \item{\code{semantic_min_sum}}{Measure of semantic distance proposed in
##'   [1]. In a nutshell, For computing semantic distance between documents A
##'   (column in x) and B (column in y), first for each term a in A the minimal
##'   distance to terms in B is computed with \code{dist_type} distance. Then,
##'   this values are summed with weights co-responding weights (\code{x}
##'   matrix). Same procedure applies to terms from B, the resulting two values
##'   are summed:
##'
##'        \deqn{DIST(A, B)=\sum_a x_{A,a}\min_b D(a,b) + \sum_b x_{B,b}\min_a D(b,a)}
##'
##'   Note that in [1] the authors weight each term by normalized IDF
##'   weight. The formulation in this package is more general. You can achieve
##'   their formula by applying "idf" \code{strans} and "l1" \code{ptrans}
##'   transformations. See examples.}
##'
##' \item{\code{semantic_min_max}}{Measure of semantic similarity proposed in
##'   [2]. The authors used the name "Relaxed Word Mover Distance" to emphasize
##'   that the measure is a lower bound of the well known "Earth Mover Distance"
##'   transportation problem. The metric is a variation of
##'   \code{semantic_min_sum} where the \code{max} is used in the last step
##'   instead of \code{sum}}
##'
##' \item{\code{adist_rwmd}}{Relaxed Word Mover Distance - same as
##'   \code{adist_semantic_min_max}.}
##' 
##' }
##'
##' @references
##'
##'   [1] Mihalcea, Rada, Courtney Corley, and Carlo Strapparava. ‘Corpus-Based
##'     and Knowledge-Based Measures of Text Semantic Similarity’. In AAAI,
##'     6:775–80, 2006.
##' 
##'   [2] Ye, Xin, Hui Shen, Xiao Ma, Razvan Bunescu, and Chang Liu. ‘From Word
##'     Embeddings to Document Similarities for Improved Information Retrieval
##'     in Software Engineering’. In Proceedings of the 38th International
##'     Conference on Software Engineering, 404–415. ICSE ’16. New York, NY,
##'     USA: ACM, 2016. doi:10.1145/2884781.2884862.
##' 
##' @param vecs Dense matrix with columns
##' @param precompute logical Weather to optimize the computation for speed and
##'     precompute individual distances. The computation is method specific bug
##'     generally should be \code{TRUE} (the default) unless memory usage is a
##'     concern.
##' @return A matrix of the distances. If \code{y=NULL}, the value is a cross
##'     distance of \code{x}.
##' @rdname aggregated_distances
##' @aliases aggr_similarities
##' @param dist_type distance to use across individual vectors in \code{vecs}
##' @inheritParams dist_cosine
##' @rdname aggregated_distances
##' @export
adist_centroid <- function(x, y = NULL, vecs, ptrans = NULL, strans = NULL,
                           by = c("primary", "secondary", "row", "column"),
                           pairwise = FALSE, precompute = !pairwise,
                           dist_type = "cosine") {
    generic_aggr_dist(x, y, vecs, ptrans, strans, by = by,
                      precompute = precompute, pairwise = pairwise, 
                      aggr_type = "CENTROID", dist_type = dist_type)
}

##' @rdname aggregated_distances
##' @export
adist_semantic_min_sum <- function(x, y = NULL, vecs, ptrans = NULL, strans = NULL,
                                   by = c("primary", "secondary", "row", "column"),
                                   pairwise = FALSE, precompute = !pairwise,
                                   dist_type = "cosine") {
    generic_aggr_dist(x, y, vecs, ptrans, strans, by = by,
                      precompute = precompute, pairwise = pairwise, 
                      aggr_type = "SEM_MIN_SUM", dist_type = dist_type)
}

##' @rdname aggregated_distances
##' @export
adist_semantic_min_max <- function(x, y = NULL, vecs, ptrans = NULL, strans = NULL,
                                   by = c("primary", "secondary", "row", "column"),
                                   pairwise = FALSE, precompute = !pairwise,
                                   dist_type = "cosine") {
    generic_aggr_dist(x, y, vecs, ptrans, strans, by = by,
                      precompute = precompute, pairwise = pairwise, 
                      aggr_type = "SEM_MIN_MAX", dist_type = dist_type)
}

## ##' @rdname aggregated_distances
## ##' @export
## adist_dtw <- function(x, y = NULL, vecs, ptrans = NULL, strans = NULL,
##                       by = c("primary", "secondary", "row", "column"),
##                       pairwise = FALSE, precompute = !pairwise,
##                       dist_type = "cosine") {
    
##     generic_aggr_dist(x, y, vecs, ptrans, strans, by = by,
##                       precompute = precompute, pairwise = pairwise, 
##                       aggr_type = "SEM_MIN_MAX", dist_type = dist_type)
## }

## tt <- matrix(runif(100), 5, 10)
## plot(dtw(tt, distance.only = F))
## plot(symmetric2)

##' @rdname aggregated_distances
##' @export
adist_rwmd <- adist_semantic_min_max

