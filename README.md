# Functional ecology and imperfect detection of species

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1065446.svg)](https://doi.org/10.5281/zenodo.1065446)

This repository contains all the materials needed to reproduce the analyses in Roth, T., Allan, E., Pearman, P. B., Amrhein, V. (2017). Functional ecology and imperfect detection of species. Accepted by *Methods in Ecology and Evolution*. These materials are presented as an R Package that contains:

- a function to simulate observed meta-community data from communities that are subject to ecological and detection filtering, 
- the analysed plant data from the Swiss Biodiversity Monitoring, 
- the values for the three analysed functional traits: specific leaf area, canopy height and seed mass for recorded species, 
- a vignette that develops the ideas behind the simulation of the meta-community, 
- a vignette that describes the workflow to estimate detection-corrected meta-community from observations using the hierarchical models implemented in the r-package umarked, and
- a vignette that describes all the analyses conducted in the manuscript.

## Installation of the package

The latest version can be downloaded directly in R using:

````
library(devtools)
install_github("TobiasRoth/detectionfilter")
library(detectionfilter)
````

## Data

The plant survey and trait data used in the manuscript are available under:

````
plantsBDM 
traitmat
````

## Vignettes

The vignette that develops the ideas behind the simulation of the meta-community that can be used to test the presented methods to account for imperfect detection on functional divesity estimates along an environmental gradient:

````
vignette("simdat", package="detectionfilter")
````

The vignette that describes the workflow to estimate a detection-corrected meta-community from observations using the hierarchical models implemented in the r-package umarked:

````
vignette("workflow", package="detectionfilter")
````

The vignette that describes all the analyses conducted in the manuscript:

````
vignette("main_analyses", package="detectionfilter")
````

## Modifications
*07-Jul-2017*: First working version of package with a first draft of all vignettes, functions and data files. [Version 0.0.1](https://github.com/TobiasRoth/detectionfilter/releases/tag/0.0.1).

*10-Aug-2017*: Edit of main_analyses.Rmd, including minor changes to graphs and clarification of text.  Some questions are inserted IN CAPS and need to be addressed. Addition of README.txt, modification of .gitignore

*11-Aug-2017*: Copy /inst/doc/main_analysis.Rmd to /vignettes/main_analyses.Rmd

*18-Aug-2017*: TOC to all vignettes. Analyses added to main_analyses.Rmd to compare functional composition and diversity measures from all observation and from observations with rare (< four observations) species removed. 

*21-Aug-2017*: Analyses added to main_analyses.Rmd to infer whether the detection effect we observe for functional packing is larger than what we expect by adding missed species at random. 

*22-Aug-2017*: Version of package submitted to Methods in Ecology and Evolution [Version 0.0.2](https://github.com/TobiasRoth/detectionfilter/releases/tag/0.0.2).

*23-Nov-2017*: Small edits to vignettes. First release of package. [Version v1.0](https://github.com/TobiasRoth/detectionfilter/releases/tag/v1.0).
