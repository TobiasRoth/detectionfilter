# Functional ecology and imperfect detection of species

Citation for submitted manuscript:

Roth, T., Allan, E., Pearman, P. B., Amrhein, V. (submitted). Functional ecology and imperfect detection of species. Submitted to Methods in Ecology and Evolution.

This repository contains all the materials needed to reproduce the analyses in Roth et al. (submitted) Functional ecology and imperfect detection of species. These materials are presented as an R Package that contains:

* a function to simulate observed meta-community data from communities that are subject to ecological and detection filtering, 
* the analysed plant data from the Swiss Biodiversity Monitoring, 
* the values for the three analysed functional traits: specific leaf area, canopy height and seed mass for recorded species, 
* a vignette that develops the ideas behind the simulation of the meta-community, 
* a vignette that describes the workflow to estimate detection-corrected meta-community from observations using the hierarchical models implemented in the r-package umarked, and
* a vignette that describes all the analyses conducted in the manuscript.

## Installation of the package

The latest version can be downloaded directly in R using:

     library(devtools)
     install_github("TobiasRoth/detectionfilter”)
     library(detectionfilter)```

## Data

The plant survey and trait data used in the manuscript are available under:

     plantsBDM 
     traitmat

##Vignettes

The vignette that develops the ideas behind the simulation of the meta-community that can be used to test the presented methods to account for imperfect detection on functional divesity estimates along an environmental gradient:
vignette("simdat", package="detectionfilter")

The vignette that describes the workflow to estimate a detection-corrected meta-community from observations using the hierarchical models implemented in the r-package umarked:

     vignette("workflow", package="detectionfilter")

The vignette that describes all the analyses conducted in the manuscript:

     vignette("main_analyses", package="detectionfilter")

## Modifications

10-Aug-2017  Edit of main_analyses.Rmd, including minor changes to graphs and clarification of text.  Some questions are inserted IN CAPS and need to be addressed. Addition of README.txt, modification of .gitignore

11-Aug-2017  Copy /inst/doc/main\_analysis.Rmd to /vignettes/main\_analyses.Rmd
