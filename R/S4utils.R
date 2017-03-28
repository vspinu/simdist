

### Checkers

check_colnames <- function(x, arg_name = deparse(substitute(x))) {
    if (is.null(colnames(x)))
        stop(arg_name, " must have valid columns names", call. = FALSE)
    colnames(x)
}

check_rownames <- function(x, arg_name = deparse(substitute(x))) {
    if (is.null(rownames(x)))
        stop(arg_name, " must have valid rownames names", call. = FALSE)
    rownames(x)
}



## Internal methods for sparse tdm and dtm manipulation

is_primary_by <- function(x, by = c("primary", "secondary", "row", "column")) {
    by <- match.arg(by)
    if (by == "row")
        primary_dim(x) == 1
    else if (by == "column")
        primary_dim(x) == 2
    else by == "primary"
}

get_norm_fun <- function(x) {
    if (!is.null(x))
        if (is.character(x))
            get(paste0("norm_", x), mode = "function")
        else if (is.function(x))
            x
        else
            stop("normalization parameter must be a character or a function")
}

.unames <- function(vec){
    if(is.factor(vec)) levels(vec)
    else {
        unames <- attr(vec, "unames")
        if (!is.null(unames))
            unames
        else
            sort(unique(vec))
    }    
}

setGeneric("primary_ix", useAsDefault = function(x) stop("Unimplimented for class ", class(x)))
setMethod("primary_ix", "matrix", function(x) 0:(ncol(x) - 1L))
setMethod("primary_ix", "TsparseMatrix", function(x) x@i)
setMethod("primary_ix", "data.frame",
          function(x) {
              vec <- primary_ids(x)
              if (is.factor(vec)) {
                  as.integer(vec) - 1L
              } else {
                  fastmatch::fmatch(vec, .unames(vec)) - 1L
              }
          })

setGeneric("secondary_ix", useAsDefault= function(x) stop("Unimplimented for class ", class(x)))
setMethod("secondary_ix", "matrix", function(x) 0:(nrow(x) - 1L))
setMethod("secondary_ix", "CsparseMatrix", function(x) x@i)
setMethod("secondary_ix", "RsparseMatrix", function(x) x@j)
setMethod("secondary_ix", "TsparseMatrix", function(x) x@j)
setMethod("secondary_ix", "data.frame",
          function(x) {
              vec <- secondary_ids(x)
              if (is.factor(vec)) {
                  as.integer(vec) - 1L
              } else {
                  fastmatch::fmatch(vec, .unames(vec)) - 1L
              }
          })

setGeneric("primary_dim", useAsDefault= function(x) stop("Unimplimented for class ", class(x)))
setMethod("primary_dim", "matrix", function(x) 2L)
setMethod("primary_dim", "CsparseMatrix", function(x) 2L)
setMethod("primary_dim", "RsparseMatrix", function(x) 1L)
setMethod("primary_dim", "TsparseMatrix", function(x) 1L)

setGeneric("secondary_dim", useAsDefault= function(x) stop("Unimplimented for class ", class(x)))
setMethod("secondary_dim", "matrix",        function(x) 2L)
setMethod("secondary_dim", "CsparseMatrix", function(x) 2L)
setMethod("secondary_dim", "RsparseMatrix", function(x) 1L)
setMethod("secondary_dim", "TsparseMatrix", function(x) 1L)

setGeneric("primary_size", useAsDefault= function(x) stop("Unimplimented for class ", class(x)))
setMethod("primary_size", "matrix", ncol)
setMethod("primary_size", "CsparseMatrix", ncol)
setMethod("primary_size", "RsparseMatrix", nrow)
setMethod("primary_size", "TsparseMatrix", nrow)
setMethod("primary_size", "data.frame", function(x) length(.unames(x[[df_var_name(x, "primary")]])))

setGeneric("secondary_size", useAsDefault= function(x) stop("Unimplimented for class ", class(x)))
setMethod("secondary_size", "matrix", nrow)
setMethod("secondary_size", "CsparseMatrix", nrow)
setMethod("secondary_size", "RsparseMatrix", ncol)
setMethod("secondary_size", "TsparseMatrix", ncol)
setMethod("secondary_size", "data.frame", function(x) length(.unames(x[[df_var_name(x, "secondary")]])))

setGeneric("primary_names", useAsDefault= function(x, ...) stop("Unimplimented for class ", class(x)))
setMethod("primary_names", "matrix", colnames)
setMethod("primary_names", "CsparseMatrix", colnames)
setMethod("primary_names", "RsparseMatrix", rownames)
setMethod("primary_names", "TsparseMatrix", rownames)
setMethod("primary_names", "data.frame", function(x) .unames(x[[df_var_name(x, "primary")]]))

setGeneric("secondary_names", useAsDefault= function(x, ...) stop("Unimplimented for class ", class(x)))
setMethod("secondary_names", "matrix", rownames)
setMethod("secondary_names", "CsparseMatrix", rownames)
setMethod("secondary_names", "RsparseMatrix", colnames)
setMethod("secondary_names", "TsparseMatrix", colnames)
setMethod("secondary_names", "data.frame", function(x) .unames(x[[df_var_name(x, "secondary")]]))


df_var_name <- function(x, type = c("primary", "secondary", "value")) {
    type <- match.arg(type)
    name <- attr(x, paste0(".", type, "_name"))
    if (!is.null(name))
        name
    else 
        names(x)[[c(primary = 1L, secondary = 2L, value = 3L)[[type]]]]
}

setGeneric("primary_ids", useAsDefault = function(x) {
    pnames <- primary_names(x)
    if (!is.null(pnames)) {
        structure(primary_ix(x) + 1L, class = "factor", levels = pnames)
    }
})
setMethod("primary_ids", "data.frame",
          function(x) x[[df_var_name(x, "primary")]])


setGeneric("secondary_ids", useAsDefault =  function(x) {
    snames <- secondary_names(x)
    if (!is.null(snames)) {
        structure(secondary_ix(x) + 1L, class = "factor", levels = snames)
    }
})
setMethod("secondary_ids", "data.frame",
          function(x) x[[df_var_name(x, "secondary")]])


setGeneric("sparse_vals", function(x) stop("Unimplimented for class ", class(x)))
setMethod("sparse_vals", "matrix", function(x) c(x))
setMethod("sparse_vals", "sparseMatrix", function(x) x@x)
setMethod("sparse_vals", "data.frame",
          function(x) x[[df_var_name(x, "value")]])


setGeneric("primary_subset", function(x, ix) stop("Unimplimented for class ", class(x)))
setMethod("primary_subset", "matrix", function(x, ix) x[, ix, drop = F])
setMethod("primary_subset", "CsparseMatrix", function(x, ix) x[, ix, drop = F])
setMethod("primary_subset", "RsparseMatrix", function(x, ix) x[ix, , drop = F])
setMethod("primary_subset", "TsparseMatrix", function(x, ix) x[ix, , drop = F])

setGeneric("secondary_subset", function(x, ix) stop("Unimplimented for class ", class(x)))
setMethod("secondary_subset", "matrix", function(x, ix) x[ix, , drop = F])
setMethod("secondary_subset", "CsparseMatrix", function(x, ix) x[ix, , drop = F])
setMethod("secondary_subset", "RsparseMatrix", function(x, ix) {
    tmp <- new("dgCMatrix",
               i = x@j, p = x@p, x = x@x,
               Dim = rev(dim(x)), Dimnames = rev(dimnames(x)))[ix,, drop = F]
    new("dgRMatrix", j = tmp@i, p = tmp@p, x = tmp@x,
        Dim = rev(dim(tmp)), Dimnames = rev(dimnames(tmp)))
})
setMethod("secondary_subset", "TsparseMatrix", function(x, ix) x[, ix, drop = F])



setGeneric("to_primary_maybe", signature = c("x"), 
           def = function(x, by) standardGeneric("to_primary_maybe"), 
           useAsDefault = function(x, by) stop("Unimplimented for class ", class(x)))
setMethod("to_primary_maybe", "matrix", function(x, by) {
    if (is_primary_by(x, by)) x else t(x)
})
setMethod("to_primary_maybe", "CsparseMatrix", function(x, by) {
    if (is_primary_by(x, by)) x else as(t(x), "RsparseMatrix")
})
setMethod("to_primary_maybe", "RsparseMatrix", function(x, by) {
    if (is_primary_by(x, by)) x else as(t(x), "CsparseMatrix")
})
setMethod("to_primary_maybe", "TsparseMatrix", function(x, by) {
    if (is_primary_by(x, by))
        x
    else 
        new("dgTMatrix", j = x@i, i = x@i, x = x@x, Dim = dim(x), Dimnames = dimnames(x))            
})
setMethod("to_primary_maybe", "data.frame", function(x, by) {
    if (is_primary_by(x, by)) x
    else {
        psv(x,
            pname = df_var_name(x, "secondary"),
            sname = df_var_name(x, "primary"),
            vname = df_var_name(x, "value"))
    }
})

