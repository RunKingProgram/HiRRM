# HiRRM -Early Access

## 1. Getting started
####	Downloading HiRRM
HiRRM can be downloaded https://github.com/RunKingProgram/HiRRM. It can be installed as a regular R package.
####	Installing HiRRM
HiRRM links to R packages Rcpp, RcppArmadillo, RcppEigen, snow, parallel,data.table, nlme and BEDMatrix. These dependencies should be installed before installing HiRRM. In addition, **HiRRM requires Hi_RRM file (Hi_RRM_linux for Linux) under your working directory**. Here is an example for installing HiRRM and all its dependencies in an R session(assuming none of the R packages other than the default has been installed):

Under Linux and MacOS
```
install.packages( c( "Rcpp", "RcppArmadillo", "RcppEigen", "snow", "parallel", "data.table", "nlme" , "BEDMatrix" ), repos = "https://cran.r-project.org/" )
system( “R CMD install HiRRM_1.0.tar.gz” )
```


## 2. Main functions
The current version of HiRRM includes two main functions:
```
coefy = Estimate_coefy(Phenotype,order) 
Hi_RRM(Genotype,coefy)
```
**Estimate_coefy** is to estimate phenotypic regression coefficients in optimal linear longitudinal trajectory with first hierarchical RRM or LS method.
**Hi_RRM** is to associate multiple phenotypic regressions with markers using EMMAX with second hierarchical mvLMM method.

#### Arguments
#### Phenotype
**Phenotype must be adjusted for non-time-dependent testing date.** Phenotype is an object class of character，which is dmu format (https://dmu.ghpc.au.dk/DMU/Doc/Current/dmuv6_guide.5.2.pdf), and it should be recorded in a text file with colname that can be read into R as a matrix. Here is an example of the header and first 9 rows for the phenotype (as "phenotye.txt" file): 
|id| trait|time|
| ---------- | :-----------:  | :-----------: |
|1 |4.9|1|
|2 |4.6|1|
|3 |8|1|
|1 |10.1|2|
|1 |14.9|3|
|2 |16.2|3|
|3 |15.2|3|
|4 |14.3|3|
...

If there is time-dependent fixed factors,such as population stratification and sex, it should be placed before anylized trait and HiRRM will adjuest it automatically
For example (as "phenotype2.txt" file):

|id|sex|trait | time|
| ---------- | :-----------:  | :-----------: | :-----------: |
|1 |1|4.9|1|
|2 |2|4.6|1|
|3 |1|8|1|
|1 |1|10.1|2|
|1 |1|14.9|3|
|2 |2|16.2|3|
|3 |1|15.2|3|
|4 |1|14.3|3|
...

#### order
An object class of numeric: the optimal growth trajectory could be chosen according to the Bayesian information criterion (BIC).
```
coefy <- Estimate_coefy(Phenotype,4) 
```

#### Genotype
An object class of character，which consists of the three files in PLINK BED format (http://www.cog-genomics.org/plink/1.9/formats#bed). For example, Genotype= “Genotype.bed”, “Genotype.bim” and “Genotype.bam” :
```
result <- Hi_RRM(“Genotype”,coefy)
```


## 3.Example
```
library(HiRRM)
library(parallel)
library(BEDMatrix)
library(data.table)
library(snow)
library(RcppArmadillo)
library(nlme)

setwd("./example")

coefy <-  Estimate_coefy("phenotype.txt",4) 
result <- Hi_RRM("Genotype",coefy)



