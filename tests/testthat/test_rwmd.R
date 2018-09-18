context("rwmd")

library(Matrix)
wv <- readRDS(system.file("test-data/wv100.rds", package = "simdist"))
tdm <- readRDS(system.file("test-data/tdm100.rds", package = "simdist"))
wv2 <- norm_l2(wv)

test_that("rwmd is computed ocrrectly on dgRMatrix", {

    tdmR <- as(t(tdm), "RsparseMatrix")
    ## !! Matrix implementation detail, subsetting Rmatrix ends into Tmatrix.
    ## class(tdmR[, 1:2])
    ## class(tdmR[1:2, ])
    tdmR10 <- as(tdmR[1:10, ], "RsparseMatrix")

    d0 <- adist_rwmd(tdm, tdm[, 1:10], wv2,  "l1", precompute = TRUE)
    d1 <- adist_rwmd(tdmR, tdmR10, wv2, "l1", precompute = TRUE)
    d2 <- adist_rwmd(tdmR, tdmR10, wv2, "l1", precompute = FALSE)
    
    expect_equal(d0, d1)
    expect_equal(d1, d2)

})

test_that("rwmd is computed correctly on dgCMatrix", {

    ## library(text2vec)
    ## rwmd <- text2vec::RWMD$new(t(wv))
    ## d0 <- rwmd$dist2(t(tdm), t(tdm[, 1:10]))
    ## saveRDS(d0, "../../inst/test-data/text2vec_rwmd100.rds")
    d0 <- readRDS(system.file("test-data/text2vec_rwmd100.rds", package = "simdist"))

    wv21 <- wv2[, 1:10]
    d1 <- adist_rwmd(tdm, tdm[, 1:10], wv2, "l1", precompute = TRUE)
    d2 <- adist_rwmd(tdm, tdm[, 1:10], wv2, "l1", precompute = FALSE)
    
    expect_equal(d1, d2)
    expect_equal(d0, unname(d1))
    expect_equal(d0, unname(d2))

})

test_that("rwmd is computed correctly on psv", {

    ## library(text2vec)
    ## rwmd <- text2vec::RWMD$new(t(wv))
    ## d0 <- rwmd$dist2(t(tdm), t(tdm[, 1:10]))
    ## saveRDS(d0, "../../inst/data/text2vec_rwmd100.rds")
    d0 <- readRDS(system.file("test-data/text2vec_rwmd100.rds", package = "simdist"))
    ## d0 <- mat2psv(d0)

    df <- mat2psv(tdm)
    df100 <- mat2psv(tdm[, 1:10])
    d1 <- adist_rwmd(df, df100, wv2, "l1", precompute = TRUE)
    d2 <- adist_rwmd(df, df100, wv2, "l1", precompute = FALSE)
    
    expect_equal(d1, d2)
    expect_equal(d0, unname(d1))
    expect_equal(d0, unname(d2))
})


test_that("rwmd works as expected with named objects", {

    set.seed(100)
    cmat1 <- simdist:::random_sparse_mat(123, 10, 25)
    cmat2 <- cmat1[rev(1:nrow(cmat1)), ]
    cmat3 <- cmat1[sample(1:nrow(cmat1)), ]

    snames <- unique(rownames(cmat1))
    vecs <- matrix(runif(50*length(snames)), 50, dimnames = list(NULL, snames))
    
    d11 <- adist_rwmd(cmat1, cmat1, vecs, "l2", dist_type = "cosine")
    d12 <- adist_rwmd(cmat1, cmat2, vecs, "l2", dist_type = "cosine")
    d13 <- adist_rwmd(cmat1, cmat3, vecs, "l2", dist_type = "cosine")
    
    expect_equal(d11, d12)
    expect_equal(d11, d13)

    d21 <- adist_rwmd(cmat2, cmat1, vecs, "l1", dist_type = "euclidean")
    d22 <- adist_rwmd(cmat2, cmat2, vecs, "l1", dist_type = "euclidean")
    d23 <- adist_rwmd(cmat2, cmat3, vecs, "l1", dist_type = "euclidean")
    expect_equal(d21, d22)
    expect_equal(d21, d23)
    
})

test_that("pairwise rwmd works as expected", {

    set.seed(200)
    cmat1 <- simdist:::random_sparse_mat(123, 10, 25)
    cmat2 <- cmat1[rev(1:nrow(cmat1)), sample(1:ncol(cmat1))]
    cmat3 <- cmat1[sample(1:nrow(cmat1)), sample(1:ncol(cmat1))]

    snames <- unique(rownames(cmat1))
    vecs <- matrix(runif(50*length(snames)), 50, dimnames = list(NULL, snames))

    d11 <- adist_rwmd(cmat1, cmat1, vecs, "l2", dist_type = "cosine")
    pd11 <- adist_rwmd(cmat1, cmat1, vecs, "l2", pairwise = T, dist_type = "cosine")
    d21 <- adist_rwmd(cmat2, cmat1, vecs, "l2", dist_type = "cosine")
    pd21 <- adist_rwmd(cmat2, cmat1, vecs, "l2", pairwise = T, dist_type = "cosine")
    d22 <- adist_rwmd(cmat2, cmat2, vecs, "l2", dist_type = "cosine")
    pd22 <- adist_rwmd(cmat2, cmat2, vecs, "l2", pairwise = T, dist_type = "cosine")
    d23 <- adist_rwmd(cmat2, cmat3, vecs, "l2", dist_type = "cosine")
    pd23 <- adist_rwmd(cmat2, cmat3, vecs, "l2", pairwise = T, dist_type = "cosine")

    expect_equal(unname(diag(d11)), pd11)
    expect_equal(unname(diag(d21)), pd21)
    expect_equal(unname(diag(d22)), pd22)
    expect_equal(unname(diag(d23)), pd23)

    d11 <- adist_rwmd(cmat1, cmat1, vecs, "l2", dist_type = "euclidean")
    pd11 <- adist_rwmd(cmat1, cmat1, vecs, "l2", pairwise = T, dist_type = "euclidean")
    d21 <- adist_rwmd(cmat2, cmat1, vecs, "l1", dist_type = "euclidean")
    pd21 <- adist_rwmd(cmat2, cmat1, vecs, "l1", pairwise = T, dist_type = "euclidean")
    d22 <- adist_rwmd(cmat2, cmat2, vecs, "l1", dist_type = "euclidean")
    pd22 <- adist_rwmd(cmat2, cmat2, vecs, "l1", pairwise = T, dist_type = "euclidean")
    d23 <- adist_rwmd(cmat2, cmat3, vecs, "l1", dist_type = "euclidean")
    pd23 <- adist_rwmd(cmat2, cmat3, vecs, "l1", pairwise = T, dist_type = "euclidean")

    expect_equal(unname(diag(d11)), pd11)
    expect_equal(unname(diag(d21)), pd21)
    expect_equal(unname(diag(d22)), pd22)
    expect_equal(unname(diag(d23)), pd23)
})

test_that("centroid works as expected with named objects", {

    set.seed(100)
    cmat1 <- simdist:::random_sparse_mat(123, 10, 25)
    cmat2 <- cmat1[rev(1:nrow(cmat1)), ]
    cmat3 <- cmat1[sample(1:nrow(cmat1)), ]

    snames <- unique(rownames(cmat1))
    vecs <- matrix(runif(50*length(snames)), 50, dimnames = list(NULL, snames))
    
    d11 <- adist_centroid(cmat1, cmat1, vecs, "l2", dist_type = "cosine")
    d12 <- adist_centroid(cmat1, cmat2, vecs, "l2", dist_type = "cosine")
    d13 <- adist_centroid(cmat1, cmat3, vecs, "l2", dist_type = "cosine")
    
    expect_equal(d11, d12)
    expect_equal(d11, d13)

    d21 <- adist_centroid(cmat2, cmat1, vecs, "l1", dist_type = "euclidean")
    d22 <- adist_centroid(cmat2, cmat2, vecs, "l1", dist_type = "euclidean")
    d23 <- adist_centroid(cmat2, cmat3, vecs, "l1", dist_type = "euclidean")
    expect_equal(d21, d22)
    expect_equal(d21, d23)
    
})

test_that("pairwise centroid works as expected", {

    set.seed(100)
    cmat1 <- simdist:::random_sparse_mat(123, 10, 25)
    cmat2 <- cmat1[rev(1:nrow(cmat1)), sample(1:ncol(cmat1))]
    cmat3 <- cmat1[sample(1:nrow(cmat1)), sample(1:ncol(cmat1))]

    snames <- unique(rownames(cmat1))
    vecs <- matrix(runif(50*length(snames)), 50, dimnames = list(NULL, snames))

    d11 <- adist_centroid(cmat1, cmat1, vecs, "l2", dist_type = "cosine")
    pd11 <- adist_centroid(cmat1, cmat1, vecs, "l2", pairwise = T, dist_type = "cosine")
    d21 <- adist_centroid(cmat2, cmat1, vecs, "l2", dist_type = "cosine")
    pd21 <- adist_centroid(cmat2, cmat1, vecs, "l2", pairwise = T, dist_type = "cosine")
    d22 <- adist_centroid(cmat2, cmat2, vecs, "l2", dist_type = "cosine")
    pd22 <- adist_centroid(cmat2, cmat2, vecs, "l2", pairwise = T, dist_type = "cosine")
    d23 <- adist_centroid(cmat2, cmat3, vecs, "l2", dist_type = "cosine")
    pd23 <- adist_centroid(cmat2, cmat3, vecs, "l2", pairwise = T, dist_type = "cosine")

    expect_equal(unname(diag(d11)), pd11)
    expect_equal(unname(diag(d21)), pd21)
    expect_equal(unname(diag(d22)), pd22)
    expect_equal(unname(diag(d23)), pd23)

    d11 <- adist_centroid(cmat1, cmat1, vecs, "l2", dist_type = "euclidean")
    pd11 <- adist_centroid(cmat1, cmat1, vecs, "l2", pairwise = T, dist_type = "euclidean")
    d21 <- adist_centroid(cmat2, cmat1, vecs, "l1", dist_type = "euclidean")
    pd21 <- adist_centroid(cmat2, cmat1, vecs, "l1", pairwise = T, dist_type = "euclidean")
    d22 <- adist_centroid(cmat2, cmat2, vecs, "l1", dist_type = "euclidean")
    pd22 <- adist_centroid(cmat2, cmat2, vecs, "l1", pairwise = T, dist_type = "euclidean")
    d23 <- adist_centroid(cmat2, cmat3, vecs, "l1", dist_type = "euclidean")
    pd23 <- adist_centroid(cmat2, cmat3, vecs, "l1", pairwise = T, dist_type = "euclidean")

    expect_equal(unname(diag(d11)), pd11)
    expect_equal(unname(diag(d21)), pd21)
    expect_equal(unname(diag(d22)), pd22)
    expect_equal(unname(diag(d23)), pd23)
})

test_that("centroid computation is correct", {
    set.seed(100)
    mat <- matrix(1:9, nrow = 3, dimnames = list(1:3, letters[1:3]))
    cmat <- as(mat, "dgCMatrix")
    vecs <- t(mat)

    ave <- cbind(a = colSums(t(vecs) * mat[, "a"]),
                 b = colSums(t(vecs) * mat[, "b"]),
                 c = colSums(t(vecs) * mat[, "c"]))/3

    cave <- as(ave, "dgCMatrix")
    cosine <- dist_cosine(cave, cave, "l2")
    euclid <- dist_euclidean(cave, cave)
    acosine <- adist_centroid(cmat, cmat, vecs, dist_type = "cosine")
    aeuclid <- adist_centroid(cmat, cmat, vecs, dist_type = "euclidean")
    expect_equal(cosine, acosine)
    expect_equal(euclid, aeuclid)    
})
