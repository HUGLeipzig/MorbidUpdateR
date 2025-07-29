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
devtools::install_github("HUGLeipzig/MorbidUpdateR")
```

## Config file

Some functions require values specified in a config.yml file. The required values are: 
- `omim_id`: your personal OMIM identifier, needed for downloading the OMIM files
- `varvis_target`: your Varvis target, needed for the URL generation, e.g. "your-university"
- `varvis_user`: the username for your Varvis API
- `varvis_password`: the password for your Varvis API (can be Unicode)
- `hgmd_csv_path`: the path and filename to your stored HGMD csv file
- `panelapp_tsv_path`: the path and filename where the PanelApp downloads should be stored
- `sysndd_tsv_path`: the path and filename where the SysNDD downloads should be stored

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

## Default CutOffs

The functions in this package have precomputed cutoffs to determine if a gene is a morbidgene or not. If _**ONE**_ of the following points is true, then a genes is considered to be a morbidgene:

- `>=4` pathogenic variants in **ClinVar**
    - column `clinvar_pathogenic_cutoff`
- `>=4` pathogenic variants in **HGMD** (= "DM" Variant)
    - column `hgmd_pathogenic_cutoff`
- `TRUE` has an **OMIM** phenotype
    - column `omim_phenotype`
- `TRUE` is green gene in **PanelApp** England _OR_ **PanelApp** Australia
    - column `panelapp`
- `TRUE` has status 'definitive' in **GenCC**
    - column `gencc`
- `TRUE` has status 'definite' in **SysNDD**
    - column `sysndd`
- `TRUE` was added **manually**
    - column `manually_added`
