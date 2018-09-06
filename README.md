# GIMLET
R package for the identification of biological modulators in context-specific gene regulation using local energy statistics

<strong>Depends:</strong>

R(>=3.3.2)

Rcpp, RcppArmadillo, infotheo

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

## Inputs of GIMLET

- **X**: a data matrix of regulators
- **Y**: a data matrix of targets
- **Z**: a data matrix of modulators
- **r**: a number of the upper and lower points of Z, that is, |U_Z| and |L_Z| in equation (4). The default value is 1.
- **delta**: a tuning parameter that indicates the proportion of neighbors in the Euclidean distance of Z. The default value is 0.3.
- **nperm**: number of permutations for calculating the p-value. The default value is 1000.
- **type**: a character string indicating what type of Y is used. One of "continuous"(default) or "binary" can be abbreviated.

## An example of synthetic data

Let denote X, Y, and Z by a regulator, its target gene, and a modulator. We show an example of synthetic data of (X,Y,Z) generated using the endogenous switching regression model where the relationship between X and Y is linear and varies with Z. For details, see the Section 3 in the manuscript.

```
library("GIMLET")

simn <- 100 # number of simulations
n <- 200 # number of samples
sim.typ <- 1 # 1:M_1, 2:M_2, 3:M_3 in Section 3
fun.typ <- 1 # 1:line, 2: quadratic, 3:cubic, 4:sine period 1/2, 5:sine period 1/8, 6:x^1/4, 7: circle, 8:step
pval_thres <- 0.05 # threshold of significance

pval <- rep(1,simn)

for(k in 1:simn){

  cat(k,"/",simn,"\n")
  set.seed(k)

  if(sim.typ==5|sim.typ==7){
    r1 <- 1/6
    r2 <- 1
  } else {
    r1 <- 0.5
    r2 <- 3
  }

  data <- sim(n=n,sim.typ=sim.typ,fun.typ=fun.typ,model=1,r1=r1,r2=r2)
  x <- data$x
  y <- data$y
  z <- data$z
  if(is.vector(x)) x <- matrix(x,ncol=1)
  if(is.vector(y)) y <- matrix(y,ncol=1)
  if(is.vector(z)) z <- matrix(z,ncol=1)

  gimlet.obj <- gimlet(x,y,z,1,0.25,1000)
  pval[k] <- gimlet.obj$global.p.value

}

mean(pval <= pval_thres)

```

## Reference
Teppei Shimamura, Yusuke Matsui, Taisuke Kajino, Satoshi Ito, Takashi Takahashi, and Satoru Miyano, GIMLET: Identifying Biological Modulators in Context-Specific Gene Regulation Using Local Energy Statistics, Accepted in proceedings of CIBB 2017. The post-proceedings version of CIBB 2017 paper is now reviewed.
