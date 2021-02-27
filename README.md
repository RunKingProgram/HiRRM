# HiRRM -Early Access

## 1. Getting started
####	Downloading HiRRM
HiRRM can be downloaded https://github.com/RunKingProgram/HiRRM. It can be installed as a regular R package.
####	Installing HiRRM
HiRRM links to R packages Rcpp, RcppArmadillo, RcppEigen, snow, parallel,data.table, nlme and BEDMatrix. These dependencies should be installed before installing HiRRM. In addition, **HiRRM requires Hi_RRM file (Hi_RRM_linux for Linux) under your working directory**. Here is an example for installing HiRRM and all its dependencies in an R session(assuming none of the R packages other than the default has been installed):


For Linux:
```
install.packages( c( "Rcpp", "RcppArmadillo", "RcppEigen", "snow", "parallel", "data.table", "nlme" , "BEDMatrix" ), repos = "https://cran.r-project.org/" )
system( “R CMD INSTALL HiRRM_1.0.tar.gz” )
```

For MacOS:
```
install.packages( c( "Rcpp", "RcppArmadillo", "RcppEigen", "snow", "parallel", "data.table", "nlme" , "BEDMatrix" ), repos = "https://cran.r-project.org/" )
system( “R CMD INSTALL HiRRM_1.0_MacOS.tar.gz” )
```
For windows OS:

We recommend Windows Subsystem for Linux (https://docs.microsoft.com/en-us/windows/wsl/about).
The installation steps are as follows:

1.Open the Microsoft Store and choose your preferred Linux distribution.
For example: Ubuntu 16.04 (https://www.microsoft.com/store/apps/9pjn388hp8c9)

2.Open the PowerShell as an administrator, enter 
```
enable windows optionalfeature - Online - featurename Microsoft Windows subsystem Linux
```
and press enter to execute the command.

3.tap following command to install R 
```
sudo apt update
sudo apt-get install r-base-core
```
4. Use "explorer.exe ." Command to management files just like Windows OS.

Now the subsystem has been installed and configured, and you can install HiRRM with the following commond:

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
|20|75031133|JAX00714218|-1.06205|-0.783078|-0.91676|0.0370952|0.036248|0.042004|31.257171|7.503e-07|
|20|78404667|JAX00181136|-1.0113|-0.736987|-0.889224|0.036135|0.035578|0.04119|29.304733|1.935e-06|
|20|78572249|JAX00181149|-0.9280|-0.6576|-0.812616|0.0330044|0.03307|0.03785|27.82017|3.96179e-06|
|20|78611877|JAX00181152|-1.1096|-0.7924|-0.97957|0.034687|0.034363|0.03954|37.4572|3.682002e-08|
|20|81285198|JAX00714779|-1.0248|-0.7564|-0.92873|0.03663|0.0354|0.041438|29.8937|1.453422e-06|

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



