# RIOM: Randomization Inference with Outcome Misclassification

## Author
Siyu Heng and Pamela A. Shaw

## Maintainer
Siyu Heng (Email: <siyuheng@sas.upenn.edu>)

## Description
**RIOM** is an R package for randomization inference with outcome misclassification, which includes the calculation of warning accuracy and sensitivity weights introduced in Heng and Shaw (2022).

Before installing this R package, please ensure that you have installed the following two R packages: **gurobi** and **mgcv**. To install package **RIOM** in R from GitHub, please run the following commands:

```
install.packages("devtools") 
library(devtools) 
install_github("siyuheng/RIOM")
```
## Reference
Heng, S. and Shaw, P. A. (2022). "A model-free and finite-population-exact framework for randomized experiments subject to outcome misclassification via integer programming." arXiv:2201.03111.
