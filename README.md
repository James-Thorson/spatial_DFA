Description
=============

SpatialDFA
* Is an R package for implementing spatial dynamic factor analysis, a model for estimating spatiotemporal density dependence in community dynamics

Background

Usage
=============

### Please contact James Thorson if you're planning to use this package

Installation Instructions
=============
This function depends on R version >=3.1.1 and a variety of other tools.

First, install the "devtools" package from CRAN

    # Install and load devtools package
    install.packages("devtools")
    library("devtools")

Second, please install the following:
* TMB (Template Model Builder): https://github.com/kaskr/adcomp
* INLA (integrated nested Laplace approximations): http://www.r-inla.org/download

Note: at the moment, TMB and INLA can be installed using the commands 

    # devtools command to get TMB from GitHub
    install_github("kaskr/adcomp/TMB") 
    # source script to get INLA from the web
    source("http://www.math.ntnu.no/inla/givemeINLA.R")  
    
Next, please install the SpatialDFA from this GitHub repository using a function in the "devtools" package:

    # Install package
    install_github("james-thorson/SpatialDFA", auth_token="918291a743d6fe53aecec3bfc4a27d2b850c688d") 
    # Load package
    library(SpatialDFA)

