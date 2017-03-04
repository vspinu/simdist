[![Build Status](https://travis-ci.org/vspinu/simdist.png?branch=master)](https://travis-ci.org/vspinu/simdist)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/simdist)](https://cran.r-project.org/package=simdist)
[![Development version](https://img.shields.io/badge/devel-0.0.1.9000-orange.svg)](https://github.com/vspinu/simdist)
[![CRAN version](http://www.r-pkg.org/badges/version/simdist)](https://cran.r-project.org/package=simdist)

== `Simdist` [R](https://www.r-project.org/) Package 

High performance distances and similarities for various dense and sparse
representations especially geared towards applications in NLP and recommender
systems.

== Supported and Planned Object Types

- `matrix` from base R
- `dgCMatrix`, `dgRMatrix` and `dgTMatrix` from Matrix package
- `simple_triplet_matrix` from `slam` package
- `data.frames` in primary-secondary-value (psv) format
- `list` of named numeric or character vectors

== Distances for 2d representations

|             | `matrix` | `dgCMatrix` | `dgRMatrix` | `dgTMatrix` | `slam` | `psv`    | `list` |
| ---:        | :---:    | :---:       | :---:-      | :---:       | :---:  | :---:    | :---:  |
| cosine      | &#10004; | &#10004;    | &#10004;    | &#10004;    |        | &#10004; |        |
| euclidean   | &#10004; | &#10004;    | &#10004;    | &#10004;    |        | &#10004; |        |
| mahalanobis |          |             |             |             |        |          |        |
| jaccard     |          |             |             |             |        |          |        |


== Distances for 3d representations

|                  | `dgCMatrix` | `dgRMatrix` | `dgTMatrix` | `slam` | `psv`    | `list` |
| ---:             | :---:       | :---:-      | :---:       | :---:  | :---:    | :---:  |
| centroid         | &#10004;    | &#10004;    | &#10004;    |        | &#10004; |        |
| semantic_min_max | &#10004;    | &#10004;    | &#10004;    |        | &#10004; |        |
| semantic_min_sum | &#10004;    | &#10004;    | &#10004;    |        | &#10004; |        |


== Transformations

`norm_l1`, `norm_l2`.






