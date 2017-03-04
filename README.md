[![Build Status](https://travis-ci.org/vspinu/simdist.png?branch=master)](https://travis-ci.org/vspinu/simdist)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/simdist)](https://cran.r-project.org/package=simdist)
[![Development version](https://img.shields.io/badge/devel-0.0.0.9000-orange.svg)](https://github.com/vspinu/simdist)
[![CRAN version](http://www.r-pkg.org/badges/version/simdist)](https://cran.r-project.org/package=simdist)

High performance distances and similarities for various dense and sparse
representations with primary focus on applications in NLP and recommender
systems.

## Supported and Planned Object Types

- `matrix` from base R
- `dgCMatrix`, `dgRMatrix` and `dgTMatrix` from Matrix package
- `simple_triplet_matrix` from `slam` package
- `data.frames` in primary-secondary-value (psv) format
- `list` of named numeric or character vectors

## Distances for 2D Representations

|               | `matrix` | `dgCMatrix` | `dgRMatrix` | `dgTMatrix` | `slam` | `psv`    | `list` |
| ---:          | :---:    | :---:       | :---:       | :---:       | :---:  | :---:    | :---:  |
| `cosine`      | &#10004; | &#10004;    | &#10004;    | &#10004;    |        | &#10004; |        |
| `euclidean`   | &#10004; | &#10004;    | &#10004;    | &#10004;    |        | &#10004; |        |
| `mahalanobis` |          |             |             |             |        |          |        |
| `jaccard`     |          |             |             |             |        |          |        |


## Aggregation Distances for 3D Representations

|                                | `dgCMatrix` | `dgRMatrix` | `dgTMatrix` | `slam` | `psv`    | `list` |
| ---:                           | :---:       | :---:       | :---:       | :---:  | :---:    | :---:  |
| `centroid`                     | &#10004;    | &#10004;    | &#10004;    |        | &#10004; |        |
| `semantic_min_max`<sup>1</sup> | &#10004;    | &#10004;    | &#10004;    |        | &#10004; |        |
| `semantic_min_sum`<sup>2</sup> | &#10004;    | &#10004;    | &#10004;    |        | &#10004; |        |

[1] More commonly known as "Relaxed Word Mover Distance" (RWMD) proposed in _Kusner et. al. [‘From Word Embeddings To Document Distances’](http://jmlr.org/proceedings/papers/v37/kusnerb15.pdf) (2015)_.


[2] Similar to RWMD measure, proposed in _Mihalcea et.al. ['Corpus-Based and Knowledge-Based Measures of Text Semantic Similarity'](https://pdfs.semanticscholar.org/1374/617e135eaa772e52c9a2e8253f49483676d6.pdf) (2006)_


## Transformations

`norm_l1`, `norm_l2`.






