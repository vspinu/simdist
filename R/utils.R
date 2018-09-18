
##' Various utilities
##'
##' \code{psv}: mark a data frame with attributes in order to track which are
##' primary, seocndary and value columns.
##' 
##' @param df data.frame
##' @param pname primary id name
##' @param sname secondary id name
##' @param vname value name
##' @param mat matrix
##' @rdname utils
##' @export
mat2psv <- function(mat) {
  smat <-
    if (inherits(mat, "sparseMatrix"))
      as(mat, "TsparseMatrix")
    else {
      dms <- dimnames(mat)
      if (is.null(dms)) dms <- list(NULL, NULL)
      new("dgTMatrix",
          i = c(row(mat)) - 1L, j = c(col(mat)) - 1L,
          x = as.double(mat),
          Dim = dim(mat), Dimnames = dms)
    }

  pri_row <- primary_dim(mat) == 1L
  pid <- if(pri_row) primary_ids(smat) else secondary_ids(smat)
  sid <- if(pri_row) secondary_ids(smat) else primary_ids(smat)
  if (is.null(pid))
    pid <- if(pri_row) primary_ix(smat) else secondary_ix(smat)
  if (is.null(sid))
    sid <- if(pri_row) secondary_ix(smat) else primary_ix(smat)
  data.frame(pid = pid, sid = sid, val = smat@x, stringsAsFactors = FALSE)
}


##' @rdname utils
##' @export
psv2mat <- function(df) {
  ## fixme: factor is slow, try fastmatch?
  pid <- factor(primary_ids(df), exclude = NULL)
  sid <- factor(secondary_ids(df), exclude = NULL)
  val <- sparse_vals(df)
  ## returns dgCMatrix, columsn are primary direction
  sparseMatrix(i = as.integer(sid), j = as.integer(pid), x = val, 
               dims = c(length(levels(sid)), length(levels(pid))),
               dimnames = list(levels(sid), levels(pid)))
}


##' @rdname utils
##' @export
psv <- function(df, pname = NULL, sname = NULL, vname = NULL) {
  stopifnot(is.data.frame(df))
  if(!is.null(pname)) {
    stopifnot(pname %in% names(df))
    attr(df, ".primary_name") <- pname
  }
  if(!is.null(sname)) {
    stopifnot(sname %in% names(df))
    attr(df, ".secondary_name") <- sname
  }
  if(!is.null(vname)) {
    stopifnot(vname %in% names(df))
    attr(df, ".value_name") <- vname
  }
  df
}



### internal utilities

stopif <- function (...) {
  n <- length(ll <- list(...))
  if (n == 0L) 
    return(invisible())
  mc <- match.call()
  for (i in 1L:n)
    if (!is.logical(r <- ll[[i]]) || anyNA(r) || any(r)) {
      ch <- deparse(mc[[i + 1]], width.cutoff = 60L)
      if (length(ch) > 1L)
        ch <- paste(ch[1L], "....")
      stop(sprintf("%s is TRUE", ch), call. = FALSE, domain = NA)
    }
  invisible()
}

random_mat <- function(seed = NULL, R = 100, C = 25){
  if(is.null(seed))
    seed <- getOption("recolab_test_seed", 17)
  set.seed(seed)
  dimnames <- list(1:R, suppressWarnings(make.names(rep_len(letters, C), unique = T)))
  matrix(sample(c(rep.int(0, 1000), 1:300), R*C, T), R, C, dimnames = dimnames)
}

random_sparse_mat <- function(seed = NULL, R = 100, C = 25){
  Matrix::Matrix(random_mat(seed, R, C), sparse = T)
}

drop_attr <- function(x, keep = c("dim", "dimnames")) {
  for (nm in names(attributes(x))) {
    if (!nm %in% keep)
      attr(x, nm) <- NULL
  }
  x
}
