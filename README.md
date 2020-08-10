# MAPCO：A mixed model for Mendelian randomization in the presence of correlated pleiotropic effects
This package provides a Mendelian randomization approach for the inference of causal relationship between traits based on summary-level GWAS data, which takes account of the horizontal pleiotropy and correlation between genetic effects on traits. 
This package provides a Mendelian randomization approach for the inference of causal relationship between traits based on summary-level GWAS data, which takes account of the horizontal pleiotropy and correlation between genetic effects on traits. It requires neither the InSIDE assumption nor the conditions on the proportion of valid instrumental variables.
## Setup
Use the following command in R to install the package:
```
library(devtools)
install_github("sqx-afk/MAPCO")
```
## Usage
```
MAPCO(by,bx,se_y,se_x,int_beta=0,int_sigma_alpha=NA,b0=0,tol=1e-8,n_iter=3000)
```
## Example 
```
library(MAPCO)
by <- data[,1]
bx <- data[,2]
se_y <- data[,3]
se_x <- data[,4]
result <- MAPCO(by,bx,se_y,se_x)
result
```
## Authors


## References

A mixed model for Mendelian randomization in the presence of correlated pleiotropic effects. 2020.(submitted)
