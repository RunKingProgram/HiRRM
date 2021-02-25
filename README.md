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
## 3.Output
The HiRRM.txt contains the results.The number of columns will depend on the number of phenotypes used for analysis. The first three columns are: chromosome numbers,base pair positions on the chromosome and snp ids. The last column contains p values from the association tests. The middle columns contain beta estimates and the variance matrix for these estimates. And QTN candidates above the Bonferroni corrected level of significance would be selected in R as matrix.

|CHR|	POS	|SNP|	BETA1	|BETA2|	BETA3|	BETA4|	BETA5|	VBETA1	|VBETA2	|VBETA3	|VBETA4|	VBETA5	|chisq|	p-value|
| ---------- | :-----------:  | :-----------: | :-----------: | ---------- | :-----------:  | :-----------: | :-----------: | ---------- | :-----------:  | :-----------: | :-----------: |
|2	|168798341|	UNC4465781|	0.87197|	0.9642221|	0.7717525	|0.9575686	|-0.04723711|	0.02376907	|0.02651553	|0.02761537	|0.03023668	|0.000999294	|41.64283	|6.96E-08|
|4	|3785255	|UNC6639515|	0.1759812|	0.8056519	|0.7291173|	0.4126044	|0.01886118	|0.02634507	|0.02924049	|0.03057526	|0.03364178	|0.00088208	|50.43038	|1.13E-09
|4	|3787560	|JAX00115433|	0.1759812|	0.8056519|	0.7291173	|0.4126044|	0.01886118|	0.02634507	|0.02924049|0.03057526	|0.03364178	|0.00088208	|50.43038	|1.13E-09|
|4	|4878592	|UNC6651698|	0.231539|	0.8650871|	0.7884397|	0.4663801	|0.01819409|	0.02694842	|0.02990495	|0.03127136	|0.03441789|	0.000882081|	52.15886	|5.00E-10|
|4	|5999840	|UNC6664472| 0.2129025|	0.8520055|	0.7457233|	0.4366588	|0.0208303	|0.02778343	|0.03082014	|0.03223382	|0.03549616|	0.000876123	|49.33987	|1.89E-09|
|4  |7033031	|UNC6677287|	0.1600859|	0.8128049|	0.6470668	|0.4120346	|0.02304759	|0.02941886	|0.03261758|	0.03412091	|0.03760379	|0.000876228|	44.15184	2.16E-08|

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

coefy <-  Estimate_coefy("phenotype.txt",4) 
result <- Hi_RRM("Genotype",coefy)



