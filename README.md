<h1>Functional ecology and imperfect detection of species</h1>

<p>Citation for submitted manuscript:</p>

<p>Roth, T., Allan, E., Pearman, P. B., Amrhein, V. (submitted). Functional ecology and imperfect detection of species. Submitted to Methods in Ecology and Evolution.</p>

<p>This repository contains all the materials needed to reproduce the analyses in Roth et al. (submitted) Functional ecology and imperfect detection of species. These materials are presented as an R Package that contains:</p>

<ul>
<li>a function to simulate observed meta-community data from communities that are subject to ecological and detection filtering, </li>
<li>the analysed plant data from the Swiss Biodiversity Monitoring, </li>
<li>the values for the three analysed functional traits: specific leaf area, canopy height and seed mass for recorded species, </li>
<li>a vignette that develops the ideas behind the simulation of the meta-community, </li>
<li>a vignette that describes the workflow to estimate detection-corrected meta-community from observations using the hierarchical models implemented in the r-package umarked, and</li>
<li>a vignette that describes all the analyses conducted in the manuscript.</li>
</ul>

<h2>Installation of the package</h2>

<p>The latest version can be downloaded directly in R using:</p>

<pre><code> library(devtools)
 install_github("TobiasRoth/detectionfilter”)
 library(detectionfilter)```
</code></pre>

<h2>Data</h2>

<p>The plant survey and trait data used in the manuscript are available under:</p>

<pre><code> plantsBDM 
 traitmat
</code></pre>

<h2>Vignettes</h2>

<p>The vignette that develops the ideas behind the simulation of the meta-community that can be used to test the presented methods to account for imperfect detection on functional divesity estimates along an environmental gradient:
vignette("simdat", package="detectionfilter")</p>

<p>The vignette that describes the workflow to estimate a detection-corrected meta-community from observations using the hierarchical models implemented in the r-package umarked:</p>

<pre><code> vignette("workflow", package="detectionfilter")
</code></pre>

<p>The vignette that describes all the analyses conducted in the manuscript:</p>

<pre><code> vignette("main_analyses", package="detectionfilter")
</code></pre>

<h2>Modifications</h2>

<p>10-Aug-2017  Edit of main_analyses.Rmd, including minor changes to graphs and clarification of text.  Some questions are inserted IN CAPS and need to be addressed. Addition of README.txt, modification of .gitignore</p>

<p>11-Aug-2017  Copy /inst/doc/main_analysis.Rmd to /vignettes/main_analyses.Rmd</p>
