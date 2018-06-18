# GIMLET
R package for the identification of biological modulators in context-specific gene regulation using local energy statistics

<strong>Depends:</strong>

R(>=3.3.2)

Rcpp, RcppArmadillo

<strong>Authors:</strong>

Teppei Shimamura

<strong>Contact:</strong>
shimamura[at]med.nagoya-u.ac.jp

## Installation:

To install GIMLET package, we start R and run the following codes on R console:

```
library(devtools)
install_github("tshimam/GIMLET")
```

## General overview

We propose a novel method, genome-wide identification of modulators using local energy statistical test (GIMLET), to identify biological modulators that contribute to transcription factor activities.
GIMLET takes data matrix (for example, expression, copy number, and mutation data matrix for n samples) of regulators, their target genes, and modulators as inputs, and calculates a new dependence measure, called local distance correlation, to compare the difference of distance correlation at low and high values of given modulators. Local distance correlation is used for modeling the relationships between regulators and their target genes at the fixed point of modulators. GIMLET calculates a statistical significance whether local distance correlation varies with modulators using a permutation-based approach.

## An example of synthetic data

Let denote X, Y, and Z by a regulator, its target gene, and a modulator. We show an example of synthetic data of (X,Y,Z) generated using the endogenous switching regression model where the relationship between X and Y is linear and varies with modulators. For details, see the Section 3 in the manuscript.

```
library("GIMLET")

n <- 200
junk <- sim1(n,type="linear")
X <- junk$x
Y <- junk$y
Z <- junk$z
obj <- gimlet(X,Y,Z)
obj$global.p.value

```

## Reference
Teppei Shimamura, Yusuke Matsui, Taisuke Kajino, Satoshi Ito, Takashi Takahashi, and Satoru Miyano, GIMLET: Identifying Biological Modulators in Context-Specific Gene Regulation Using Local Energy Statistics, Accepted in CIBB 2017. The post-proceedings version of CIBB 2017 paper will be published in Lecture Notes in Bioinformatics.
