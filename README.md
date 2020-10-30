# HiRRM -Early Access

Package: HiRRM<br>
Type: Package<br>
Title: Hi-RRM<br>
Version: 0.9<br>
Date: 2020-10-08<br>
Author: Yuxin Song<br>
Maintainer: Runqing Yang <runqingyang@cafs.ac.cn><br>
Description: Fast Genome-wide Association Analysis for Dynamic Traits Using Hierachical Mixed-Model<br>
License: GPL (>= 2)<br>
Imports: Rcpp (>= 1.0.5)<br>
LinkingTo: Rcpp,RcppArmadillo,RcppEigen,snow,parallel,data.table,nlme,BEDMatrix<br>
Built: R 4.0.0; x86_64-apple-darwin19.4.0; 2020-10-08 08:29:35 UTC; unix<br>

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
