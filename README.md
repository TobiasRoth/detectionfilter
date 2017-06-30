# Functional ecology and imperfect detection of species

Citation for submitted manuscript:

Roth, T., Allan, E., Pearman, P. B., Amrhein, V. (submitted). Functional ecology and imperfect detection of species. Submitted to Methods in Ecology and Evolution.

This repository contains all the materials needed to reproduce Roth *et al.* (submitted) Functional ecology and imperfect detection of species.  These materials are presented as an R Package which contains a function to simulate observed meta-community data from communities that are subject to ecological and detection filtering, the analysed plants data from the Swiss Biodiversity Monitoring, the values for the three functional traits specific leaf area, canopy height and seed mass for recorded species, a vignette that develops the ideas behind the simulation of the meta-community, a vignette that describes the workflow to estimate detection-corrected meta-community from observations using the hierarchical models implemented in the r-package umarked, and a vignette that describes all the analyses conducted in the manuscript.



## Installation of the package
The latest version can be downloaded directly in R using:
```
library(devtools)
install_github("TobiasRoth/detectionfilter")
library(detectionfilter)
```

## Data
The plants survey and trait data used in this manuscript are available under:
```
plantsBDM 
traitmat
```

## Vignettes
The vignette that develops the ideas behind the simulation of the meta-community that can be used to test the presented methods to account for imperfect detection on functional divesity estimates along an environmental gradient:
```
vignette("simdat", package="detectionfilter")
```

The vignette that describes the workflow to estimate a detection-corrected meta-community from observations using the hierarchical models implemented in the r-package umarked: 
```
vignette("workflow", package="detectionfilter")
```

The vignette that describes all the analyses conducted in the manuscript:
```
vignette("main_analyses", package="detectionfilter")
```
