# HiRRM -Early Access

## 1. Getting started
####	Downloading HiRRM
HiRRM can be downloaded https://github.com/YuxinSong-prog/HiRRM. It can be installed as a regular R package.
####	Installing HiRRM
HiRRM links to R packages Rcpp, RcppArmadillo, RcppEigen, snow, data.table, nlme and BEDMatrix. These dependencies should be installed before installing HiRRM. In addition, HiRRM requires Hi_RRM file (https://github.com/YuxinSong-prog/HiRRM) under your run directory. Here is an example for installing HiRRM and all its dependencies in an R session(assuming none of the R packages other than the default has been installed):
```
install.packages(c("Rcpp", "RcppArmadillo", "RcppEigen", "snow", "data.table", "nlme" , "BEDMatrix"), repos = "https://cran.r-project.org/")
system(“R CMD install HiRRM_1.0.tgz”)
```
## 2. Main functions
At release HiRRM will include two main functions:
```
coefy = Fit_Curve(y,inputorders) 
Hi_RRM(plinkfilename,coefy)
```
#### Arguments
#### y
Phenotype should be either saved as a data frame in R, or recorded in a text file that can be read into R as a data frame. The rows of the data frame represent different individuals, the columns represent different time, and the columns represent variable.
#### inputorders
Order of polynomial you need to fit.
#### plinkfilename
An object class of character: the filename before the suffix of plink BED files. The three plink file must have a same filename, for example, “plinkfilename.bed”, “plinkfilename.bim” and “plinkfilename.bam”.

## 3.Example
```
library(HiRRM)
library(BEDMatrix)
library(data.table)
library(snow)
library(nlme)

path1<-"/Users/songyuxin/Desktop/CFWmice_binary"
setwd(path1)
y=as.matrix(fread("mvpheno"))
coefy = Fit_Curve(y,3) 
result = Hi_RRM("fill_cfw",coefy)
```
