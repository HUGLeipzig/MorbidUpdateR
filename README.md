---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->



# MorbidUpdateR

<!-- badges: start -->
<!-- badges: end -->

The goal of MorbidUpdateR is to conveniently update our Morbid Genes Panel. The user 
can simply run the functions without editing the variables, as the default values 
already point to the correct and relevant files and paths. 

For more flexibility and for testing purposes, the user can also specify new paths, 
versions and files (which is, however, not recommended, as the alternative input 
values have not been extensively tested yet).

So in summary, just run the functions every month for an efficient update of the 
Morbid Genes Panel. 

## Installation

You can install the development version from [GitHub](https://github.com/) with:


```r
# install.packages("devtools")
#devtools::install_github("RJauss/Morbid-UpdateR")
```

## Example

This is a basic example which shows you how to solve a common problem:


```r
#library(MorbidUpdateR)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so:


```r
#summary(cars)
```

You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. You could also use GitHub Actions to re-render `README.Rmd` every time you push. An example workflow can be found here: <https://github.com/r-lib/actions/tree/master/examples>.

You can also embed plots, for example:



In that case, don't forget to commit and push the resulting figure files, so they display on GitHub and CRAN.
