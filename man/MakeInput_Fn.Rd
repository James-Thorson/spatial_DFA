\name{MakeInput_Fn}
\alias{MakeInput_Fn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to bundle inputs to TMB based on data-frame and either default or advanced settings
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MakeInput_Fn(Version, Nfactors, DF, inla_spde, Kappa_Type = "Constant", ObsModel = NULL, Include_Omega = TRUE, Include_Epsilon = TRUE, EncounterFunction = 2, Include_Phi = TRUE, Include_Rho = TRUE, Use_REML = FALSE, X_ik = NULL, YearSet = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Version}{
%%     ~~Describe \code{Version} here~~
}
  \item{Nfactors}{
%%     ~~Describe \code{Nfactors} here~~
}
  \item{DF}{
%%     ~~Describe \code{DF} here~~
}
  \item{inla_spde}{
%%     ~~Describe \code{inla_spde} here~~
}
  \item{Kappa_Type}{
%%     ~~Describe \code{Kappa_Type} here~~
}
  \item{ObsModel}{
%%     ~~Describe \code{ObsModel} here~~
}
  \item{Include_Omega}{
%%     ~~Describe \code{Include_Omega} here~~
}
  \item{Include_Epsilon}{
%%     ~~Describe \code{Include_Epsilon} here~~
}
  \item{EncounterFunction}{
%%     ~~Describe \code{EncounterFunction} here~~
}
  \item{Include_Phi}{
%%     ~~Describe \code{Include_Phi} here~~
}
  \item{Include_Rho}{
%%     ~~Describe \code{Include_Rho} here~~
}
  \item{Use_REML}{
%%     ~~Describe \code{Use_REML} here~~
}
  \item{X_ik}{
%%     ~~Describe \code{X_ik} here~~
}
  \item{YearSet}{
%%     ~~Describe \code{YearSet} here~~
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
function (Version, Nfactors, DF, inla_spde, Kappa_Type = "Constant", 
    ObsModel = NULL, Include_Omega = TRUE, Include_Epsilon = TRUE, 
    EncounterFunction = 2, Include_Phi = TRUE, Include_Rho = TRUE, 
    Use_REML = FALSE, X_ik = NULL, YearSet = NULL) 
{
    if (is.null(YearSet)) 
        YearSet = min(DF[, "year"]):max(DF[, "year"])
    if (is.null(ObsModel)) {
        ObsModel = ifelse(all(is.integer(ObsModel[, "catch"])), 
            0, 1)
    }
    Nyears = length(YearSet)
    Nsites = length(unique(DF[, "sitenum"]))
    Nspecies = length(levels(DF[, "spp"]))
    Nknots = mesh$n
    Nobs = nrow(DF)
    Nfactors_input = ifelse(Nfactors == 0, 1, Nfactors)
    Options_vec = c(ObsModel = ObsModel, Include_Omega = Include_Omega, 
        Include_Epsilon = Include_Epsilon, EncounterFunction = EncounterFunction)
    if (is.null(X_ik)) {
        X_ik = matrix(0, nrow = Nobs, ncol = Nspecies)
        X_ik[cbind(1:Nobs, as.numeric(DF[, "spp"]))] = 1
    }
    else {
        if (nrow(X_ik) != Nobs) 
            stop("Check design matrix input")
    }
    if (Version == "spatial_dfa_v3") 
        TmbData = list(n_obs = Nobs, n_sites = Nsites, n_years = Nyears, 
            n_knots = Nknots, n_species = Nspecies, n_factors = Nfactors_input, 
            n_cov = ncol(X_ik), c_i = DF[, "catch"], p_i = as.numeric(DF[, 
                "spp"]) - 1, s_i = DF[, "sitenum"] - 1, t_i = match(DF[, 
                "year"], YearSet) - 1, X_ik = X_ik, G0 = inla_spde$param.inla$M0, 
            G1 = inla_spde$param.inla$M1, G2 = inla_spde$param.inla$M2)
    if (Version == "spatial_dfa_v4") 
        TmbData = list(Options_vec = Options_vec, n_obs = Nobs, 
            n_sites = Nsites, n_years = Nyears, n_knots = Nknots, 
            n_species = Nspecies, n_factors = Nfactors_input, 
            n_cov = ncol(X_ik), c_i = DF[, "catch"], p_i = as.numeric(DF[, 
                "spp"]) - 1, s_i = DF[, "sitenum"] - 1, t_i = match(DF[, 
                "year"], YearSet) - 1, X_ik = X_ik, G0 = inla_spde$param.inla$M0, 
            G1 = inla_spde$param.inla$M1, G2 = inla_spde$param.inla$M2)
    if (Version \%in\% c("spatial_dfa_v10", "spatial_dfa_v9", "spatial_dfa_v8b", 
        "spatial_dfa_v8", "spatial_dfa_v7", "spatial_dfa_v6", 
        "spatial_dfa_v5")) 
        TmbData = list(Options_vec = Options_vec, n_obs = Nobs, 
            n_sites = Nsites, n_years = Nyears, n_knots = Nknots, 
            n_species = Nspecies, n_factors = Nfactors_input, 
            n_cov = ncol(X_ik), c_i = DF[, "catch"], p_i = as.numeric(DF[, 
                "spp"]) - 1, s_i = DF[, "sitenum"] - 1, t_i = match(DF[, 
                "year"], YearSet) - 1, X_ik = X_ik, Rotation_jj = diag(Nfactors_input), 
            G0 = inla_spde$param.inla$M0, G1 = inla_spde$param.inla$M1, 
            G2 = inla_spde$param.inla$M2)
    if (Version == "spatial_dfa_v2") 
        TmbParams = list(logkappa_j = rep(log(100), TmbData$n_factors), 
            alpha_j = rep(0, TmbData$n_factors), phi_j = rep(0, 
                TmbData$n_factors), loglambda_j = rep(log(1), 
                TmbData$n_factors), rho_j = rep(0.2, TmbData$n_factors), 
            Psi_val = rnorm(TmbData$n_factors * TmbData$n_species - 
                TmbData$n_factors * (TmbData$n_factors - 1)/2), 
            gamma_k = rep(1, TmbData$n_cov), log_sigma_p = rep(log(1), 
                TmbData$n_species), Epsilon_input = array(0, 
                dim = c(TmbData$n_knots, TmbData$n_factors, TmbData$n_years)), 
            Omega_input = array(0, dim = c(TmbData$n_knots, TmbData$n_factors)), 
            delta_i = rep(0, TmbData$n_obs))
    if (Version == "spatial_dfa_v3") 
        TmbParams = list(logkappa_j = rep(log(100), TmbData$n_factors), 
            alpha_j = rep(0, TmbData$n_factors), phi_j = rep(0, 
                TmbData$n_factors), loglambda_j = rep(log(1), 
                TmbData$n_factors), rho_j = rep(0.2, TmbData$n_factors), 
            L_val = rnorm(TmbData$n_factors * TmbData$n_species - 
                TmbData$n_factors * (TmbData$n_factors - 1)/2), 
            gamma_k = rep(1, TmbData$n_cov), log_sigma_p = rep(log(1), 
                TmbData$n_species), Epsilon_input = array(0, 
                dim = c(TmbData$n_knots, TmbData$n_factors, TmbData$n_years)), 
            Omega_input = array(0, dim = c(TmbData$n_knots, TmbData$n_factors)), 
            delta_i = rep(0, TmbData$n_obs))
    if (Version \%in\% c("spatial_dfa_v7", "spatial_dfa_v6", "spatial_dfa_v5", 
        "spatial_dfa_v4")) 
        TmbParams = list(logkappa_j = rep(log(100), TmbData$n_factors), 
            alpha_j = rep(0, TmbData$n_factors), phi_j = rep(0, 
                TmbData$n_factors), loglambda_j = rep(log(1), 
                TmbData$n_factors), rho_j = rep(0.2, TmbData$n_factors), 
            L_val = rnorm(TmbData$n_factors * TmbData$n_species - 
                TmbData$n_factors * (TmbData$n_factors - 1)/2), 
            gamma_k = rep(1, TmbData$n_cov), log_sigma_p = rep(log(1), 
                TmbData$n_species), zinfl_pz = matrix(0, nrow = TmbData$n_species, 
                ncol = 2), Epsilon_input = array(0, dim = c(TmbData$n_knots, 
                TmbData$n_factors, TmbData$n_years)), Omega_input = array(0, 
                dim = c(TmbData$n_knots, TmbData$n_factors)), 
            delta_i = rep(0, TmbData$n_obs))
    if (Version \%in\% c("spatial_dfa_v10", "spatial_dfa_v9", "spatial_dfa_v8b", 
        "spatial_dfa_v8")) 
        TmbParams = list(logkappa_jz = array(log(100), dim = c(TmbData$n_factors, 
            2)), alpha_j = rep(0, TmbData$n_factors), phi_j = rep(0, 
            TmbData$n_factors), loglambda_j = rep(log(1), TmbData$n_factors), 
            rho_j = rep(0.2, TmbData$n_factors), L_val = rnorm(TmbData$n_factors * 
                TmbData$n_species - TmbData$n_factors * (TmbData$n_factors - 
                1)/2), gamma_k = rep(1, TmbData$n_cov), log_sigma_p = rep(log(1), 
                TmbData$n_species), zinfl_pz = matrix(0, nrow = TmbData$n_species, 
                ncol = 2), Epsilon_input = array(0, dim = c(TmbData$n_knots, 
                TmbData$n_factors, TmbData$n_years)), Omega_input = array(0, 
                dim = c(TmbData$n_knots, TmbData$n_factors)), 
            delta_i = rep(0, TmbData$n_obs))
    Random = c("Omega_input", "Epsilon_input", "delta_i")
    if (Use_REML == TRUE) 
        Random = c(Random, "gamma_k")
    Map = list()
    Map[["alpha_j"]] = factor(rep(NA, length(TmbParams[["alpha_j"]])))
    Map[["rho_j"]] = factor(rep(1, length(TmbParams[["rho_j"]])))
    Map[["loglambda_j"]] = factor(rep(1, length(TmbParams[["loglambda_j"]])))
    if (Version \%in\% c("spatial_dfa_v7", "spatial_dfa_v6", "spatial_dfa_v5", 
        "spatial_dfa_v4")) {
        if (Kappa_Type == "Constant") 
            Map[["logkappa_j"]] = factor(rep(1, length(TmbParams[["logkappa_j"]])))
        if (Kappa_Type == "Omega_vs_Epsilon") 
            stop("Not implemented")
    }
    if (Version \%in\% c("spatial_dfa_v10", "spatial_dfa_v9", "spatial_dfa_v8b", 
        "spatial_dfa_v8")) {
        if (Kappa_Type == "Constant") 
            Map[["logkappa_jz"]] = factor(array(1, dim = dim(TmbParams[["logkappa_jz"]])))
        if (Kappa_Type == "Omega_vs_Epsilon") 
            Map[["logkappa_jz"]] = factor(outer(rep(1, TmbData$n_factors), 
                c(1, 2)))
    }
    if (Options_vec["ObsModel"] == 0) {
        Map[["zinfl_pz"]] = factor(rep(NA, prod(dim(TmbParams[["zinfl_pz"]]))))
    }
    if (Options_vec["ObsModel"] == 1) {
        Map[["delta_i"]] = factor(rep(NA, length(TmbParams[["delta_i"]])))
    }
    if (Include_Phi == FALSE | Nfactors == 0) {
        Map[["phi_j"]] = factor(rep(NA, length(TmbParams[["phi_j"]])))
        TmbParams[["phi_j"]] = rep(0, length(TmbParams[["phi_j"]]))
    }
    if (Include_Rho == FALSE | Nfactors == 0) {
        Map[["rho_j"]] = factor(rep(NA, length(TmbParams[["rho_j"]])))
        TmbParams[["rho_j"]] = rep(0, length(TmbParams[["rho_j"]]))
    }
    if (Options_vec["Include_Omega"] == FALSE | Nfactors == 0) {
        Map[["Omega_input"]] = factor(array(NA, dim = dim(TmbParams[["Omega_input"]])))
        TmbParams[["Omega_input"]] = array(0, dim = dim(TmbParams[["Omega_input"]]))
        Map[["loglambda_j"]] = factor(rep(NA, length(TmbParams[["loglambda_j"]])))
        TmbParams[["loglambda_j"]] = rep(25, length(TmbParams[["loglambda_j"]]))
    }
    if (Options_vec["Include_Epsilon"] == FALSE | Nfactors == 
        0) {
        Map[["Epsilon_input"]] = factor(array(NA, dim = dim(TmbParams[["Epsilon_input"]])))
        TmbParams[["Epsilon_input"]] = array(0, dim = dim(TmbParams[["Epsilon_input"]]))
        Map[["loglambda_j"]] = factor(rep(NA, length(TmbParams[["loglambda_j"]])))
        TmbParams[["loglambda_j"]] = rep(25, length(TmbParams[["loglambda_j"]]))
    }
    if ((Options_vec["Include_Omega"] == FALSE & Options_vec["Include_Epsilon"] == 
        FALSE) | Nfactors == 0) {
        TmbData$Options_vec[c("Include_Omega", "Include_Epsilon")] = 0
        Map[["logkappa_jz"]] = factor(array(NA, dim = dim(TmbParams[["logkappa_jz"]])))
        Map[["rho_j"]] = factor(rep(NA, length(TmbParams[["rho_j"]])))
        Map[["L_val"]] = factor(rep(NA, length(TmbParams[["L_val"]])))
    }
    if (!is.na(Options_vec["EncounterFunction"]) && Options_vec["EncounterFunction"] == 
        2) {
        Map[["zinfl_pz"]] = factor(cbind(1:TmbData$n_species, 
            rep(NA, TmbData$n_species)))
    }
    Return = list(TmbData = TmbData, TmbParams = TmbParams, Random = Random, 
        Map = Map)
    return(Return)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
