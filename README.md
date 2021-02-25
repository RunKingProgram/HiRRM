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
Best_curve=Fit_Curve(Phenotype,maxorder,fixeffect=NULL)
coefy = Estimate_coefy(Phenotype,inputorders,fixeffect=NULL) 
Hi_RRM(plinkfilename,coefy)
```
**Fit_Curve** functions is used to fit the population means of longitudinal phenotypes

**Estimate_coefy** functions is used to estimate phenotypic regression coefficients

**Hi_RRM** function is second hierarchical mvLMM, associate multiple phenotypic regressions with markers using EMMAX. 


#### Arguments
#### Phenotype
Phenotype should be either saved as a matrix in R, or recorded in a text file that can be read into R as a matrix. **Phenotype must be adjusted for non-time-dependent testing date.** Here is an example of the header and first 9 rows for the phenotype (as "phenotye.txt" file): 
|id| time| trait|
| ---------- | :-----------:  | :-----------: |
|1 |1|4.9|
|2 |1|4.6|
|3 |1|8|
|1 |2|10.1|
|1 |3|14.9|
|2 |3|16.2|
|3 |3|15.2|
|4 |3|14.3|
...
```
Phenotype <- read.table("phenotype.txt",head=T)
Best_curve=Fit_Curve(Phenotype,maxorder,fixeffect=NULL)
```
The results of **Fit_Curve** function are 
```
_order_0    R^2: 0 ; AIC: 122.6429 ; BIC: 123.4387 
_order_1    R^2: 0.97146 ; AIC: 85.522 ; BIC: 86.71569 
_order_2    R^2: 0.971475 ; AIC: 87.51621 ; BIC: 89.1078 
_order_3    R^2: 0.9999998 ; AIC: -40.5251 ; BIC: -38.53562 
_order_4    R^2: 0.9999999 ; AIC: -44.48917 ; BIC: -42.1018 
_order_5    R^2: 0.9999999 ; AIC: -43.36619 ; BIC: -40.58092 
_order_6    R^2: 0.9999999 ; AIC: -42.24417 ; BIC: -39.06101 
_order_7    R^2: 0.9999999 ; AIC: -40.85443 ; BIC: -37.27337 
_order_8    R^2: 0.9999999 ; AIC: -39.2315 ; BIC: -35.25255 
_order_9    R^2: 0.9999999 ; AIC: -39.87008 ; BIC: -35.49323  
_order_9 r2 is max  
_order_4 BIC is min  
_order_4 AIC is min  
```
and regression coefficients of the best growth curves are save as a list in R:
```
$r2
 Le[, 1:which(r2all == max(r2all))]1  Le[, 1:which(r2all == max(r2all))]2 
                        88.719767426                         88.226011856 
 Le[, 1:which(r2all == max(r2all))]3  Le[, 1:which(r2all == max(r2all))]4 
                         0.402920926                        -18.914937832 
 Le[, 1:which(r2all == max(r2all))]5  Le[, 1:which(r2all == max(r2all))]6 
                        -0.036471601                          0.006312249 
 Le[, 1:which(r2all == max(r2all))]7  Le[, 1:which(r2all == max(r2all))]8 
                        -0.017551812                         -0.001878643 
 Le[, 1:which(r2all == max(r2all))]9 Le[, 1:which(r2all == max(r2all))]10 
                         0.016968062                         -0.070140721 

$BIC
Le[, 1:which(BICall == min(BICall))]1 Le[, 1:which(BICall == min(BICall))]2 
                          88.72016734                           88.20360783 
Le[, 1:which(BICall == min(BICall))]3 Le[, 1:which(BICall == min(BICall))]4 
                           0.40424517                          -18.95483978 
Le[, 1:which(BICall == min(BICall))]5 
                          -0.03676707 

$AIC
Le[, 1:which(AICall == min(AICall))]1 Le[, 1:which(AICall == min(AICall))]2 
                          88.72016734                           88.20360783 
Le[, 1:which(AICall == min(AICall))]3 Le[, 1:which(AICall == min(AICall))]4 
                           0.40424517                          -18.95483978 
Le[, 1:which(AICall == min(AICall))]5 
                          -0.03676707 

```

#### maxorder
The max order of Legendre polynomial used to fit the population means of longitudinal phenotypes, default by 9.
#### fixeffect
fixeffect is the column number of time-dependent fixed factors,such as population stratification and sex, default by NULL.
For example (as "phenotype2.txt" file):

|id|sex| time| trait|
| ---------- | :-----------:  | :-----------: | :-----------: |
|1 |1|1|4.9|
|2 |2|1|4.6|
|3 |1|1|8|
|1 |1|2|10.1|
|1 |1|3|14.9|
|2 |2|3|16.2|
|3 |1|3|15.2|
|4 |1|3|14.3|
...

```
fixeffect = 2
Best_curve = Fit_Curve(Phenotype,fixeffect)
inputorders = 4
coefy = Estimate_coefy(Phenotype,inputorders,fixeffect) 
```

#### inputorders
An object class of numeric: the order of polynomial you need to fit.

We chose the optimal growth trajectory according to the Bayesian information criterion (BIC).
```
inputorders = 4
coefy = Estimate_coefy(Phenotype,inputorders) 
```

#### plinkfilename
An object class of character: the filename of PLINK BED files. The three PLINK files must have a same filename, for example: “plinkfilename.bed”, “plinkfilename.bim” and “plinkfilename.bam”.
For example:
```
result <- Hi_RRM("geno",coefy)
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
Phenotype <- read.table("phenotype.txt",head=T)
Best_curve <- Fit_Curve(Phenotype)


coefy <-  Fit_Curve(Phenotype,4) 
result <- Hi_RRM("geno",coefy)



