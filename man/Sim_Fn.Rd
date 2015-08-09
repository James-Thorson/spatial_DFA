\name{Sim_Fn}
\alias{Sim_Fn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to simulate data sets for testing SDFA model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Sim_Fn(n_species, n_years, n_stations = 20, phi = NULL, n_factors = 2, SpatialScale = 0.1, SD_O = 0.5, SD_E = 0.2, SD_extra = 0.1, rho = 0.8, logMeanDens = 1, Lmat = NULL, Loc = NULL, RandomSeed = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n_species}{
%%     ~~Describe \code{n_species} here~~
}
  \item{n_years}{
%%     ~~Describe \code{n_years} here~~
}
  \item{n_stations}{
%%     ~~Describe \code{n_stations} here~~
}
  \item{phi}{
%%     ~~Describe \code{phi} here~~
}
  \item{n_factors}{
%%     ~~Describe \code{n_factors} here~~
}
  \item{SpatialScale}{
%%     ~~Describe \code{SpatialScale} here~~
}
  \item{SD_O}{
%%     ~~Describe \code{SD_O} here~~
}
  \item{SD_E}{
%%     ~~Describe \code{SD_E} here~~
}
  \item{SD_extra}{
%%     ~~Describe \code{SD_extra} here~~
}
  \item{rho}{
%%     ~~Describe \code{rho} here~~
}
  \item{logMeanDens}{
%%     ~~Describe \code{logMeanDens} here~~
}
  \item{Lmat}{
%%     ~~Describe \code{Lmat} here~~
}
  \item{Loc}{
%%     ~~Describe \code{Loc} here~~
}
  \item{RandomSeed}{
%%     ~~Describe \code{RandomSeed} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (n_species, n_years, n_stations = 20, phi = NULL, n_factors = 2, 
    SpatialScale = 0.1, SD_O = 0.5, SD_E = 0.2, SD_extra = 0.1, 
    rho = 0.8, logMeanDens = 1, Lmat = NULL, Loc = NULL, RandomSeed = NA) 
{
    if (!is.na(RandomSeed)) 
        set.seed(RandomSeed)
    if (is.null(Lmat)) {
        Lmat = matrix(rnorm(n_factors * n_species), nrow = n_species, 
            ncol = n_factors)
        for (i in 1:ncol(Lmat)) {
            Lmat[seq(from = 1, to = i - 1, length = i - 1), i] = 0
            if (Lmat[, i][which.max(abs(Lmat[, i]))] < 0) {
                Lmat[, i] = -1 * Lmat[, i]
            }
        }
    }
    if (is.null(phi)) 
        phi = rnorm(n_factors, mean = 0, sd = 1)
    Beta = rep(logMeanDens, n_species)
    if (is.null(Loc)) 
        Loc = cbind(x = runif(n_stations, min = 0, max = 1), 
            y = runif(n_stations, min = 0, max = 1))
    model_O <- RMgauss(var = SD_O^2, scale = SpatialScale)
    model_E <- RMgauss(var = SD_E^2, scale = SpatialScale)
    Omega = matrix(NA, nrow = n_stations, ncol = n_factors)
    for (i in 1:n_factors) {
        Omega[, i] = RFsimulate(model = model_O, x = Loc[, "x"], 
            y = Loc[, "y"])@data[, 1]
    }
    Epsilon = array(NA, dim = c(n_stations, n_factors, n_years))
    for (i in 1:n_factors) {
        Epsilon[, i, 1] = RFsimulate(model = model_E, x = Loc[, 
            "x"], y = Loc[, "y"])@data[, 1]
        for (t in 2:n_years) {
            Epsilon[, i, t] = rho * Epsilon[, i, t - 1] + RFsimulate(model = model_E, 
                x = Loc[, "x"], y = Loc[, "y"])@data[, 1]
        }
    }
    Psi = array(NA, dim = c(n_stations, n_factors, n_years))
    for (i in 1:n_factors) {
        for (t in 1:n_years) {
            Psi[, i, t] = phi[i] * rho^t + Epsilon[, i, t] + 
                Omega[, i]/(1 - rho)
        }
    }
    Theta = array(NA, dim = c(n_stations, n_species, n_years))
    for (s in 1:n_stations) {
        for (t in 1:n_years) {
            Theta[s, , t] = Lmat \%*\% Psi[s, , t]
        }
    }
    DF = NULL
    for (s in 1:n_stations) {
        for (p in 1:n_species) {
            for (t in 1:n_years) {
                Tmp = c(sitenum = s, spp = p, year = t, catch = rpois(1, 
                  lambda = exp(Theta[s, p, t] + logMeanDens + 
                    SD_extra * rnorm(1))), waterTmpC = 0)
                DF = rbind(DF, Tmp)
            }
        }
    }
    DF = data.frame(DF, row.names = NULL)
    DF[, "spp"] = factor(letters[DF[, "spp"]])
    if (n_species > 26) 
        stop("problem with using letters")
    Sim_List = list(DF = DF, Psi = Psi, Lmat = Lmat, phi = phi, 
        Loc = Loc, Omega = Omega, Epsilon = Epsilon, Theta = Theta, 
        Psi = Psi)
    return(Sim_List)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
