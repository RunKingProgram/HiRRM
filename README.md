# HiRRM -Early Access

## 1. Getting started
####	Downloading HiRRM
HiRRM can be downloaded https://github.com/RunKingProgram/HiRRM. It can be installed as a regular R package.
####	Installing HiRRM
HiRRM links to R packages Rcpp, RcppArmadillo, RcppEigen, snow, parallel,data.table, nlme and BEDMatrix. These dependencies should be installed before installing HiRRM. In addition, **HiRRM requires Hi_RRM file (Hi_RRM_linux for Linux) under your working directory**. Here is an example for installing HiRRM and all its dependencies in an R session(assuming none of the R packages other than the default has been installed):



```
install.packages( c( "Rcpp", "RcppArmadillo", "RcppEigen", "snow", "parallel", "data.table", "nlme" , "BEDMatrix" ), repos = "https://cran.r-project.org/" )
system( “R CMD INSTALL HiRRM_1.0.tar.gz” )
```


## 2. Usage

```
coefy <- Estimate_coefy(Phenotype,Order) 
result <- Hi_RRM(Genotype,coefy)
```
**Estimate_coefy** is to estimate phenotypic regression coefficients in individual dynamic trajectory with the first hierarchical RRM or separate LS methods;
**Hi_RRM** is to associate multiple phenotypic regressions with markers using  an efficient mvLMM association analysis.

#### Arguments
#### Phenotype
An object class of character. Phenotype file contains variables ID, Covariates such as Population stratefication and Sex etc., Penotype values and Age in turn.
For example :
|ID|Sex|Phenotype| Age|
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

#### Order
An object class of numeric or formula，corrently,   which is taken at the order of optimal Legendre polynomial used to fit  population dynamic trajectory.


#### Genotype
An object class of character，which consists of three PLINK BED files with the same name. For example, Genotype.bed, Genotype.bim and Genotype.fam.


## 3.Output files
There are the two outputs from the HiRRM:  Phenotypic regressions（Regs.phenotype）and Association tests(allasso.txt and sigasso.txt).
For example, sigasso.txt:
|Chr|Pos|SNP|Beta 1|Beta 2|Beta 3|Vb 1|Vb 2|Vb 3|Chisq|P-value|
| ---------- | :-----------:  | :-----------: | :-----------: | ---------- | :-----------:  |:-----------: | ---------- | :-----------:  |:-----------: | ---------- | 
|20|75031133|JAX00714218|-1.06205415750014|-0.783078892482977|-0.91674820636336|0.0370341928042952|0.0362580842230648|0.0420037812357834|31.2571863373371|7.5038614691156e-07|
|20|78404667|JAX00181136|-1.01134786297868|-0.736987320440311|-0.889211741117024|0.0361902398083435|0.0355768512663628|0.0411202987071879|29.3047509700833|1.93236776746685e-06|
|20|78572249|JAX00181149|-0.928019777343459|-0.657628090809899|-0.812655620988216|0.0330065398293644|0.0330525657163367|0.0378037420837085|27.8201799643207|3.96177177489989e-06|
|20|78611877|JAX00181152|-1.10963838844999|-0.792408730060862|-0.97953174951917|0.0346848296947177|0.0343638986067993|0.0395464758809984|37.4574233774212|3.68209502507202e-08|
|20|81285198|JAX00714779|-1.0248959067511|-0.756450445779237|-0.928577770914673|0.0364703697603463|0.0358474434647184|0.0414384133519158|29.8930451663197|1.45342849152292e-06|

## 4.Example
```
library(HiRRM)
library(parallel)
library(BEDMatrix)
library(data.table)
library(snow)
library(RcppArmadillo)
library(nlme)

setwd("./example")
Phenotype <- "phenotype.txt"
Genotype <- "Genotype"
Order <- 2

coefy <-  Estimate_coefy(Phenotype,Order) 
result <- Hi_RRM(Genotype,coefy)



