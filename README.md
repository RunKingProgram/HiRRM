# HiRRM -Early Access

## 1. Getting started
####	Downloading HiRRM
HiRRM can be downloaded https://github.com/RunKingProgram/HiRRM. It can be installed as a regular R package.
####	Installing HiRRM
HiRRM links to R packages Rcpp, RcppArmadillo, RcppEigen, snow, parallel,data.table, nlme and BEDMatrix. These dependencies should be installed before installing HiRRM. In addition, **HiRRM requires Hi_RRM file (downloaded from https://github.com/YuxinSong-prog/HiRRM) under your working directory**. Here is an example for installing HiRRM and all its dependencies in an R session(assuming none of the R packages other than the default has been installed):

Under Linux 
```
install.packages( c( "Rcpp", "RcppArmadillo", "RcppEigen", "snow", "parallel", "data.table", "nlme" , "BEDMatrix" ), repos = "https://cran.r-project.org/" )
system( “R CMD install HiRRM.tar.gz” )
```
Under MacOS
```
install.packages( c( "Rcpp", "RcppArmadillo", "RcppEigen", "snow", "parallel", "data.table", "nlme" , "BEDMatrix" ), repos = "https://cran.r-project.org/" )
system( “R CMD install HiRRM_1.0.tgz” )
```

## 2. Main functions
The current version of HiRRM includes two main functions:
```
coefy = Fit_Curve(y,inputorders) 
Hi_RRM(plinkfilename,coefy)
```
#### Arguments
#### y
Phenotype should be either saved as a matrix in R, or recorded in a text file that can be read into R as a matrix. Here is an example of the header and first 9 rows for the phenotype: 

|time| trait| id|
| ---------- | :-----------:  | :-----------: |
|1 |4.9| 1|
|1 |4.6| 2|
|1 |8| 3|
|2 |10.1| 1|
|3 |14.9| 1|
|3 |16.2| 2|
|3 |15.2| 3|
|3 |14.3| 4|
...

#### inputorders
An object class of numeric: the order of polynomial you need to fit.
#### plinkfilename
An object class of character: the filename of PLINK BED files. The three PLINK files must have a same filename, for example: “plinkfilename.bed”, “plinkfilename.bim” and “plinkfilename.bam”.

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
y <- read.table("phenotype.txt",head=T)
coefy <-  Fit_Curve(y,4) 
result <- Hi_RRM("geno",coefy)
```
