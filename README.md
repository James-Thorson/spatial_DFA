Description
=============

SpatialDFA
* An R package for implementing spatial dynamic factor analysis, a model for estimating spatiotemporal density dependence in community dynamics.

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
    install_github("james-thorson/Spatial_DFA") 
    # Load package
    library(Spatial_DFA)

Works to read and cite if using the software
=============
* Thorson, J.T., Ianelli, J.N., Larsen, E., Ries, L., Scheuerell, M.D., Szuwalski, C. & Zipkin, E.F. (2016) Joint dynamic species distribution models: a tool for community ordination and spatiotemporal monitoring. _Global Ecology and Biogeography_. 25: 1144-1158. Abastract available [here](http://onlinelibrary.wiley.com/doi/10.1111/geb.12464/abstract).

Further information
=============
* Thorson, J.T., Fonner, R., Haltuch, M., Ono, K. & Winker, H. (2016) Accounting for spatiotemporal variation and fisher targeting when estimating abundance from multispecies fishery data. _Canadian Journal of Fisheries and Aquatic Sciences_. 73: 1-14
* Thorson, J.T., Scheuerell, M.D., Shelton, A.O., See, K.E., Skaug, H.J. & Kristensen, K. (2015) Spatial factor analysis: a new tool for estimating joint species distributions and correlations in species range. _Methods in Ecology and Evolution_, 6: 627–637.
* Thorson, J.T., Skaug, H.J., Kristensen, K., Shelton, A.O., Ward, E.J., Harms, J.H. & Benante, J.A. (2014) The importance of spatial models for estimating the strength of density dependence. _Ecology_, 96: 1202–1212.


