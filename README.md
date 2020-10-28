# HiRRM -Early Access
Package: HiRRM
Type: Package
Title: Fast Genome-wide Association Analysis for Dynamic Traits Using
        Hierachical Mixed-Model
Version: 1.0
Date: 2020-10-08
Author: Yuxin Song
Maintainer: Runqing Yang <runqingyang@cafs.ac.cn>
Description: Fast Genome-wide Association Analysis for Dynamic Traits Using Hierachical Mixed-Model
License: GPL (>= 2)
Imports: Rcpp (>= 1.0.5)
LinkingTo: Rcpp,RcppArmadillo,RcppEigen,snow,parallel,data.table,nlme,BEDMatrix
Built: R 4.0.0; x86_64-apple-darwin19.4.0; 2020-10-08 08:29:35 UTC; unix

At release HiRRM will include two main functions:
coefy = Fit_Curve(y,inputorders) and Hi_RRM(plinkfilename,coefy)
