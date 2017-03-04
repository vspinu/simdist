context("distance2d")

test_that("distances work as expected", {

    cmat <- simdist:::random_sparse_mat(1, 100, 200)
    dmat <- as.matrix(cmat)

    d1 <- as.matrix(proxy::dist(dmat, method = "cosine", by_rows = F))
    d2 <- dist_cosine(cmat, ptrans = "l2")
    d3 <- dist_cosine(dmat, ptrans = "l2")
    d4 <- dist_cosine(unname(cmat),, "l2")
    d5 <- dist_cosine(unname(dmat),, "l2")
    expect_equal(d1, d2)
    expect_equal(d1, d3)
    expect_equal(unname(d1), d4)
    expect_equal(unname(d1), d5)

    d1 <- as.matrix(proxy::dist(dmat, method = "euclidean", by_rows = F))
    d2 <- dist_euclidean(cmat)
    d3 <- dist_euclidean(dmat)
    d4 <- dist_euclidean(unname(dmat))
    d5 <- dist_euclidean(unname(cmat))
    expect_equal(d1, d2)
    expect_equal(d1, d3)
    expect_equal(unname(d1), d4)
    expect_equal(unname(d1), d5)
    
})

test_that("matrix distances work as expected", {

    dmat1 <- simdist:::random_mat(12, 10, 25)
    dmat2 <- dmat1[rev(1:nrow(dmat1)), ]
    dmat3 <- dmat1[sample(1:nrow(dmat1)), ]
    udmat <- dmat1
    rownames(udmat) <- NULL

    cmat <- as(dmat1, "dgCMatrix")
    dp <- simdist:::drop_attr(proxy::dist(dmat1, dmat1, method = "cosine", by_rows = F))
    dp["j", "j"] <- 1
    dc <- dist_cosine(cmat, , "l2")
    d1 <- dist_cosine(dmat1, dmat1, "l2")
    d2 <- dist_cosine(dmat1, dmat2, "l2")
    d3 <- dist_cosine(dmat1, dmat3, "l2")
    expect_equal(d1, dp)
    expect_equal(d1, dc)
    expect_equal(d1, d2)
    expect_equal(d1, d3)

    d4 <- dist_cosine(unname(dmat1), unname(dmat1), "l2")
    d5 <- dist_cosine(udmat, udmat, "l2")
    expect_equal(unname(d1), d4)
    expect_equal(d1, d5)
    
    cmat <- as(dmat1, "dgCMatrix")
    dc <- dist_euclidean(cmat, cmat, "l2")
    d1 <- dist_euclidean(dmat1, dmat1, "l2")
    d2 <- dist_euclidean(dmat1, dmat2, "l2")
    d3 <- dist_euclidean(dmat1, dmat3, "l2")
    expect_equal(d1, dc)
    expect_equal(d1, d2)
    expect_equal(d1, d3)

    d4 <- dist_euclidean(unname(dmat1), unname(dmat1), "l2")
    d5 <- dist_euclidean(udmat, udmat, "l2")
    expect_equal(unname(d1), d4)
    expect_equal(d1, d5)

    dp <- simdist:::drop_attr(proxy::dist(dmat1, dmat1, method = "euclidean", by_rows = F))
    d1 <- dist_euclidean(dmat1, dmat2)
    d2 <- dist_euclidean(dmat1, dmat3)
    expect_equal(d1, dp)
    expect_equal(d1, d2)

})

test_that("Csparse distances work as expected", {

    set.seed(100)
    cmat1 <- simdist:::random_sparse_mat(12, 10, 25)
    cmat2 <- cmat1[rev(1:nrow(cmat1)), ]
    cmat3 <- cmat1[sample(1:nrow(cmat1)), ]
    ucmat <- cmat1
    rownames(ucmat) <- NULL

    dmat1 <- as.matrix(cmat1)
    dp <- simdist:::drop_attr(proxy::dist(dmat1, dmat1, method = "cosine", by_rows = F))
    dp["j", "j"] <- 1
    d1 <- dist_cosine(cmat1, cmat2, "l2")
    d2 <- dist_cosine(cmat1, cmat3, "l2")
    expect_equal(dp, d1)
    expect_equal(dp, d2)
    d4 <- dist_cosine(unname(cmat1), unname(cmat1), "l2")
    d5 <- dist_cosine(ucmat, ucmat, "l2")
    expect_equal(unname(d1), d4)
    expect_equal(d1, d5)

    d1 <- simdist:::drop_attr(proxy::dist(dmat1, dmat1, method = "euclidean", by_rows = F))
    d2 <- dist_euclidean(cmat1, cmat2)
    d3 <- dist_euclidean(cmat1, cmat3)
    expect_equal(d1, d2)
    expect_equal(d1, d3)
    d4 <- dist_euclidean(unname(cmat1), unname(cmat1))
    d5 <- dist_euclidean(ucmat, ucmat)
    expect_equal(unname(d1), d4)
    expect_equal(d1, d5)
    
})

test_that("Tsparse distances work as expected", {

    library(Matrix)
    tmat1 <- as(t(simdist:::random_sparse_mat(123, 10, 25)), "TsparseMatrix")
    tmat2 <- tmat1[, rev(1:ncol(tmat1))]
    tmat3 <- tmat1[, sample(1:ncol(tmat1))]
    dmat1 <- as.matrix(tmat1)
    utmat <- tmat1
    colnames(utmat) <- NULL

    d1 <- simdist:::drop_attr(proxy::dist(dmat1, dmat1, method = "cosine", by_rows = T))
    d1[c("h", "q"), c("h", "q")] <- 1
    d11 <- dist_cosine(tmat1, tmat1, "l2")
    d2 <- dist_cosine(tmat1, tmat2, "l2")
    d3 <- dist_cosine(tmat1, tmat3, "l2")
    expect_equal(d1, d11)
    expect_equal(d1, d2)
    expect_equal(d1, d3)
    d4 <- dist_cosine(unname(tmat1), unname(tmat1), "l2")
    d5 <- dist_cosine(utmat, utmat, "l2")
    expect_equal(unname(d1), d4)
    expect_equal(d1, d5)

    d1 <- simdist:::drop_attr(proxy::dist(dmat1, dmat1, method = "euclidean", by_rows = T))
    d11 <- dist_euclidean(tmat1, tmat1)
    d2 <- dist_euclidean(tmat1, tmat2)
    d3 <- dist_euclidean(tmat1, tmat3)
    expect_equal(d1, d11)
    expect_equal(d1, d2)
    expect_equal(d1, d3)    
    d4 <- dist_euclidean(unname(tmat1), unname(tmat1))
    d5 <- dist_euclidean(utmat, utmat)
    expect_equal(unname(d1), d4)
    expect_equal(d1, d5)

})

test_that("TripleDF distances work as expected with named objects", {

    library(Matrix)
    cmat <- simdist:::random_sparse_mat(1234, 10, 2)
    dmat <- as.matrix(cmat)

    psv1 <- mat2psv(cmat)
    psv2 <- psv1[rev(1:nrow(psv1)), ]
    psv3 <- psv1[sample(1:ncol(psv1)), ]

    d1 <- simdist:::drop_attr(proxy::dist(dmat, dmat, method = "cosine", by_rows = F))
    d11 <- dist_cosine(psv1, psv1, "l2")
    d2 <- dist_cosine(psv1, psv2, "l2")
    d3 <- dist_cosine(psv1, psv3, "l2")
    expect_equal(c(d1), d11$dist)
    expect_equal(c(d1), d2$dist)
    expect_equal(c(d1), d3$dist)

    d1 <- simdist:::drop_attr(proxy::dist(dmat, dmat, method = "euclidean", by_rows = F))
    d11 <- dist_euclidean(psv1, psv1)
    d2 <- dist_euclidean(psv1, psv2)
    d3 <- dist_euclidean(psv1, psv3)
    expect_equal(c(d1), d11$dist)
    expect_equal(c(d1), d2$dist)
    expect_equal(c(d1), d3$dist)    

})
