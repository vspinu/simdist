
##' Transformations across primary and secondary dimensions
##'
##' @param x object 
##' @param by name of the dimension
##' @param na.rm logical
##' @name transform
NULL

.inf2zero <- function(x) {
    x[is.infinite(x)] <- 0
    x
}

## there is not direct dgCMatrix -> dgRMatrix coercion in Matrix; waf?
.coercion_class <- function (x) {
    switch(class(x),
           dgRMatrix = "RsparseMatrix",
           class(x))
}

.norm <- function(x, ptrans, strans) {
    x2 <- x
    if (!is.null(ptrans)) {
        fun <- get_norm_fun(ptrans)
        x2 <- fun(x, by = "primary")
    }
    if (!is.null(strans)) {
        fun <- get_norm_fun(strans)
        x2 <- fun(x, by = "secondary")
    }
    if (!identical(class(x), class(x2)))
        x2 <- as(x2, .coercion_class(x))
    x2
}

.df_norm <- function(x, by, fun) {
    primary <- is_primary_by(x, by)
    id <- if (primary) primary_ids(x) else secondary_ids(x)
    norm <- .inf2zero(ave(sparse_vals(x), id, FUN = fun))
    x[[df_var_name(x, "value")]] <- sparse_vals(x) * norm
    x
}

.mat_norm <- function(m, by, row_fun, col_fun) {
    primary <- is_primary_by(m, by)
    by_row <- if (primary) primary_dim(m) == 1 else secondary_dim(m) == 1

    norm_vec <-
        if (by_row) row_fun(m)
        else col_fun(m)

    if (by_row) {
        if (inherits(m, "CsparseMatrix"))
            ## this doesn't work for RsparseMatrix :(
            Diagonal(x = norm_vec) %*% m
        else if(inherits(m, "RsparseMatrix")) {
            m@x <- m@x * rep(norm_vec, diff(m@p))
            m
        } else
            m * norm_vec
    } else {
        if (inherits(m, "sparseMatrix"))
            m %*% Diagonal(x = norm_vec)
        else
            m * rep(norm_vec, each = nrow(m))
    }
}

##' @rdname transform
##' @export
norm_l1 <- function(x, by = c("primary", "secondary", "row", "column"), na.rm = FALSE) {
    if (is.data.frame(x))
        .df_norm(x, by, function(x) 1/sum(abs(x), na.rm = na.rm))
    else 
        .mat_norm(x, by,
                  function(x) .inf2zero(1/rowSums(abs(x), na.rm = na.rm)),
                  function(x) .inf2zero(1/colSums(abs(x), na.rm = na.rm)))
}

##' @rdname transform
##' @export
norm_l2 <- function(x, by = c("primary", "secondary", "row", "column"), na.rm = FALSE) {
    if (is.data.frame(x))
        .df_norm(x, by, function(x) 1/sqrt(sum(x^2, na.rm = na.rm)))
    else 
        .mat_norm(x, by,
                  function(x) .inf2zero(1/sqrt(rowSums(x^2, na.rm = na.rm))), 
                  function(x) .inf2zero(1/sqrt(colSums(x^2, na.rm = na.rm))))
}
