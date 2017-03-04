context("transformers")

library(Matrix)

test_that("l1 and l2 normalization works for various object classes", {

    cmat <- random_sparse_mat(1234)
    dmat <- as.matrix(cmat)
    rmat <- as(t(dmat), "RsparseMatrix")

    cl1 <- norm_l1(cmat)
    dl1 <- norm_l1(dmat)
    rl1 <- norm_l1(rmat)

    expect_equal(dl1, as.matrix(cl1))
    expect_equal(dl1, as.matrix(t(rl1)))

    cl2 <- norm_l2(cmat)
    dl2 <- norm_l2(dmat)
    rl2 <- norm_l2(rmat)

    expect_equal(dl2, as.matrix(cl2))
    expect_equal(dl2, as.matrix(t(rl2)))

})

test_that("l1 and l2 normalization returns object of input class", {

    cmat <- random_sparse_mat(1234)
    dmat <- as.matrix(cmat)
    rmat <- as(t(dmat), "RsparseMatrix")

    dl1 <- norm_l1(dmat)
    cl1 <- norm_l1(cmat)
    rl1 <- norm_l1(rmat)

    expect_equal(class(dl1), class(dmat))
    expect_equal(class(cl1), class(cmat))
    expect_equal(class(rl1), class(rmat))

})

