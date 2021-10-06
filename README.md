---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->



# MorbidUpdateR


The goal of MorbidUpdateR is to conveniently update our Morbid Genes Panel. The user 
can simply run the functions without editing the variables, as the default values 
already point to the correct and relevant files and paths. 

For more flexibility and for testing purposes, the user can also specify new paths, 
versions and files (which is, however, not recommended, as the alternative input 
values have not been extensively tested yet).

So in summary, just run the functions every month for an efficient update of the 
Morbid Genes Panel. 

## Installation

You can install the development version from [GitHub](https://github.com/HUGLeipzig/MorbidUpdateR) with:


```r
# install.packages("devtools")
GitHubToken = "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"
devtools::install_github("HUGLeipzig/MorbidUpdateR", auth_token = GitHubToken)
#> Skipping install of 'MorbidUpdateR' from a github remote, the SHA1 (dee4f84f) has not changed since last install.
#>   Use `force = TRUE` to force installation
```

## Getting started

The package contains several function, most of them can be run "as is" without 
specifying variables and paths:


```r
library(MorbidUpdateR)

# generate directory and download all the relevant files 
StartNewVersion()

# add/edit the downloaded files automatically
Add_all()

# build and save the panel
Build_MorbidGenesPanel()
```

That's basically it. For a more comprehensive list of variables and how to manipulate 
your input/output files, just run `?function` in the command line


