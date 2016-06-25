
#' Build data input for spatial dynamic factor analysis (SDFA)
#'
#' \code{MakeInput_Fn} builds a tagged list of inputs for TMB
#'
#' @param Version a version number (see example for current default).
#' @param Nfactors The number of dynamic factors used to approximate spatio-temporal variation
#' @param DF a data frame of data where each row is a unique sample and with the following columns
#' \describe{
#'   \item{catch}{the observation for each sample}
#'   \item{year}{a column of years for each sample}
#'   \item{spp}{a factor specifying the species for that sample}
#'   \item{sitenum}{a column of years for each sample}
#'   \item{PredTF}{a vector of 0s or 1st, stating whether a sample is used when fitting the model (PredTF=0) or when evaluating its predictive performance (PredTF=1)}
#'   \item{TowID}{a vector coding the unit of overdispersion, e.g., tow or vessel (OPTIONAL)}
#' }
#' @param loc_xy Locations for each station
#' @param method The method for approximating spatial variation (options: "grid" or "mesh")
#' @param Nobsfactors The number of factors used to approximate overdispersion among levels of \code{TowID} (Default=0)
#' @param Kappa_Type Whether the decorrelation distance is constant for all factors or not (Default="Constant")
#' @param ObsModel The observation model used
#' @param Include_Omega Whether to estimate purely spatial variation (Default=TRUE)
#' @param Include_Epsilon, Whether to incldue spatio-temporal variation (Default=TRUE)
#' @param EncounterFunction, The link between density and encounter probability; 0=two-parameter logistic relationship; 1=two-parameter logistic relationship with less-than-one saturation; 2=one-parameter saturating relationship (Default=1)
#' @param Correlated_Overdispersion, Whether to estimate overdispersion (only possible if TowID is present in \code{DF}
#' @param Include_Phi Whether to estimate each factor in equilibrium (FALSE) or with a fixed offset from equilibrium (TRUE), Default=TRUE
#' @param Include_Rho Whether to estimate the magnitude of temporal autocorrelation (Default=TRUE)
#' @param Use_REML Whether to use REML estimation
#' @param X_ik A matrix specifying measured variables that affect catchability for each sample (Default is turned off)
#' @param X_nl A matrix specifying measured variables that affect density for each sample (Default is an intercept)
#' @param X_ntl An array specifying measured variables that affect density and vary over time (Default is off)
#' @param a_n A vector giving the area associated with mesh vertex used when calculating indices of abundance (Default is even weighting of sample locations)
#' @param YearSet A vector of levels of the \code{year} input that should be modeled (Default is to infer from the year column of DF)
#' @param IndependentTF Whether the spatio-temporal variation (IndependentTF[1]) or overdispersion (IndependentTF[2]) is independent among species (Default is both are correlated)
#' @param CheckForBugs Whether to check inputs for obvious problems (Default=TRUE)
#' @param CorrGroup_pp A matrix filled with integers, for post hoc testing of differences in among-species correlations

#' @return Tagged list containing inputs to function \code{Build_TMB_Fn()}
#' \describe{
#'   \item{TmbData}{A tagged list of data inputs}
#'   \item{TmbParams}{A tagged list with input parameters}
#'   \item{Random}{A character vector specifying which parameters are treated as random effects}
#'   \item{Map}{A parameter map, used to turn off or mirror parameter values}
#' }

#' @export
MakeInput_Fn = function( Version, Nfactors, DF, loc_xy, method="mesh", Nobsfactors=0, Kappa_Type="Constant", ObsModel=NULL, Aniso=FALSE, Include_Omega=TRUE, Include_Epsilon=TRUE, EncounterFunction=2, Correlated_Overdispersion=FALSE, Include_Phi=TRUE, Include_Rho=TRUE, Use_REML=FALSE, X_ik=NULL, X_nl=NULL, X_ntl=NULL, a_n=NULL, YearSet=NULL, IndependentTF=c(FALSE,FALSE), CheckForBugs=TRUE, CorrGroup_pp=NULL, ...){
                                                                   
  # Calculate spde inputs
  if( require(INLA)==FALSE ) stop("Must install INLA from: source('http://www.math.ntnu.no/inla/givemeINLA.R')")
  # Build SPDE object using INLA
  inla_mesh = inla.mesh.create( loc_xy )  # loc_samp  ;  ,max.edge.data=0.08,max.edge.extra=0.2
  inla_spde = INLA::inla.spde2.matern(inla_mesh, alpha=2)

  # 2D AR1 grid
  dist_grid = dist(loc_xy, diag=TRUE, upper=TRUE)
  grid_size_km = sort(unique(dist_grid))[1]
  M0 = as( ifelse(as.matrix(dist_grid)==0, 1, 0), "dgTMatrix" )
  M1 = as( ifelse(as.matrix(dist_grid)==grid_size_km, 1, 0), "dgTMatrix" )
  M2 = as( ifelse(as.matrix(dist_grid)==sqrt(2)*grid_size_km, 1, 0), "dgTMatrix" )
  grid_list = list("M0"=M0, "M1"=M1, "M2"=M2, "grid_size_km"=grid_size_km)

  # Infer default values for inputs
  if( method=="mesh" ) Nknots = inla_mesh$n
  if( method=="grid" ) Nknots = nrow(loc_xy)
  if( is.null(YearSet) ) YearSet = min(DF[,'year']):max(DF[,'year'])
  if( is.null(ObsModel) ){
    ObsModel = ifelse( all(is.integer(DF[,'catch'])), 0, 1 )
  }
  if( is.null(a_n) ) a_n = rep(0,Nknots)
  # Species Grouping matrix
  if( is.null(CorrGroup_pp) ){
    CorrGroup_pp=matrix(0,nrow=length(levels(DF[,'spp'])),ncol=length(levels(DF[,'spp'])))
  }else{
    if( any(CorrGroup_pp!=0 & CorrGroup_pp!=1) ) stop("CorrGroup_pp can only contain 0 and 1")
    if( any(dim(CorrGroup_pp)!=length(levels(DF[,'spp']))) ) stop("CorrGroup_pp must be a square matrix with dimension Nspecies")
  }

  # Check for inconsistent inputs
  if( method=="grid" & Aniso==TRUE ){
    Aniso = FALSE
    message("Switching Aniso=FALSE because the 2D AR1 grid is isotropic")
  }

  # Options_vec
  Options_vec=c( "ObsModel"=ObsModel, "Include_Omega"=Include_Omega, "Include_Epsilon"=Include_Epsilon, "EncounterFunction"=EncounterFunction, "Correlated_Overdispersion"=ifelse(Nobsfactors==0,0,1), "AnisoTF"=Aniso, "Method"=switch(method,"mesh"=0,"grid"=1) )

  # Data size
  Nyears = length(YearSet)
  Nsites = length(unique(DF[,'sitenum']))
  Nspecies = length(unique(DF[,'spp']))
  Nobs = nrow(DF)
  Nfactors_input = ifelse( Nfactors==0, 1, Nfactors )
  Nobsfactors_input = ifelse( Nobsfactors==0, 1, Nobsfactors )
  if( Options_vec["Correlated_Overdispersion"]==1 ){
    if( !("TowID" %in% names(DF)) ) stop("with correlated observations, TowID must be a column in DF")
    Nsamples = length(unique(DF[,'TowID']))
  }else{
    if( !("TowID" %in% names(DF)) ) DF = cbind(DF, "TowID"=1)
    Nsamples = 1
  }
  if( !("PredTF_i" %in% names(DF)) ){
    DF = cbind(DF, "PredTF_i"=0)
  }else{
    if( !(all(DF[,'PredTF_i']==0 | DF[,'PredTF_i']==1)) ) stop("PredTF_i must be either 0 or 1")
  }

  # By default, catchability design matrix is turned off
  if( is.null(X_ik) ){
    X_ik = matrix(0, nrow=Nobs, ncol=1)
  }else{
    if( nrow(X_ik)!=Nobs ) stop("Check catchability design matrix input: X_ik")
  }

  # By default, spatial design matrix is an intercept for each species
  # This design matrix is not used in V14-V17
  if( is.null(X_nl) ){
    X_nl = matrix(1, nrow=Nknots, ncol=1)
  }else{
    if( !(Version%in%c("spatial_dfa_v13")) ) stop("X_nl is only used in version spatial_dfa_v13")
    if( nrow(X_nl)!=Nknots ) stop("Check spatial matrix input: X_sp")    
    if( all(sapply(X_nl,MARGIN=2,FUN=var)>0) & CheckForBugs==TRUE ) stop("You almost certainly want to add an intercept to the density matrix X_nl")
  }

  # By default, spatio-temporal design matrix is an intercept for each species
  # This design matrix is not used in V1-V13
  if( is.null(X_ntl) ){
    X_ntl = array(1, dim=c(Nknots,Nyears,1))
  }else{
    if( !(Version%in%c("spatial_dfa_v18","spatial_dfa_v17","spatial_dfa_v16","spatial_dfa_v15","spatial_dfa_v14")) ) stop("X_ntl is only used in version spatial_dfa_v14 and higher")
    if( dim(X_ntl)[1] != Nknots ) stop("Check spatial matrix input: X_ntl")    
  }

  # Data
  if(Version=="spatial_dfa_v3") TmbData = list("n_obs"=Nobs, "n_sites"=Nsites, "n_years"=Nyears, "n_knots"=Nknots, "n_species"=Nspecies, "n_factors"=Nfactors_input, "n_cov"=ncol(X_ik), "c_i"=DF[,'catch'], "p_i"=as.numeric(DF[,'spp'])-1, "s_i"=DF[,'sitenum']-1, "t_i"=match(DF[,'year'],YearSet)-1, "X_ik"=X_ik, "G0"=inla_spde$param.inla$M0, "G1"=inla_spde$param.inla$M1, "G2"=inla_spde$param.inla$M2 )
  if(Version=="spatial_dfa_v4") TmbData = list("Options_vec"=Options_vec, "n_obs"=Nobs, "n_sites"=Nsites, "n_years"=Nyears, "n_knots"=Nknots, "n_species"=Nspecies, "n_factors"=Nfactors_input, "n_cov"=ncol(X_ik), "c_i"=DF[,'catch'], "p_i"=as.numeric(DF[,'spp'])-1, "s_i"=DF[,'sitenum']-1, "t_i"=match(DF[,'year'],YearSet)-1, "X_ik"=X_ik, "G0"=inla_spde$param.inla$M0, "G1"=inla_spde$param.inla$M1, "G2"=inla_spde$param.inla$M2 )
  if(Version%in%c("spatial_dfa_v10","spatial_dfa_v9","spatial_dfa_v8b","spatial_dfa_v8","spatial_dfa_v7","spatial_dfa_v6","spatial_dfa_v5")) TmbData = list("Options_vec"=Options_vec, "n_obs"=Nobs, "n_sites"=Nsites, "n_years"=Nyears, "n_knots"=Nknots, "n_species"=Nspecies, "n_factors"=Nfactors_input, "n_cov"=ncol(X_ik), "c_i"=DF[,'catch'], "p_i"=as.numeric(DF[,'spp'])-1, "s_i"=DF[,'sitenum']-1, "t_i"=match(DF[,'year'],YearSet)-1, "X_ik"=X_ik, "Rotation_jj"=diag(Nfactors_input), "G0"=inla_spde$param.inla$M0, "G1"=inla_spde$param.inla$M1, "G2"=inla_spde$param.inla$M2 )
  if(Version%in%c("spatial_dfa_v11")) TmbData = list("Options_vec"=Options_vec, "n_obs"=Nobs, "n_samples"=Nsamples, "n_obsfactors"=Nobsfactors_input, "n_sites"=Nsites, "n_years"=Nyears, "n_knots"=Nknots, "n_species"=Nspecies, "n_factors"=Nfactors_input, "n_cov"=ncol(X_ik), "c_i"=DF[,'catch'], "m_i"=as.numeric(DF[,'TowID'])-1, "p_i"=as.numeric(DF[,'spp'])-1, "s_i"=DF[,'sitenum']-1, "t_i"=match(DF[,'year'],YearSet)-1, "X_ik"=X_ik, "Rotation_jj"=diag(Nfactors_input), "G0"=inla_spde$param.inla$M0, "G1"=inla_spde$param.inla$M1, "G2"=inla_spde$param.inla$M2 )
  if(Version%in%c("spatial_dfa_v12")) TmbData = list("Options_vec"=Options_vec, "n_obs"=Nobs, "n_samples"=Nsamples, "n_obsfactors"=Nobsfactors_input, "n_sites"=Nsites, "n_years"=Nyears, "n_knots"=Nknots, "n_species"=Nspecies, "n_factors"=Nfactors_input, "n_cov"=ncol(X_ik), "c_i"=DF[,'catch'], "m_i"=as.numeric(DF[,'TowID'])-1, "p_i"=as.numeric(DF[,'spp'])-1, "s_i"=DF[,'sitenum']-1, "t_i"=match(DF[,'year'],YearSet)-1, "X_ik"=X_ik, "a_n"=a_n, "Rotation_jj"=diag(Nfactors_input), "G0"=inla_spde$param.inla$M0, "G1"=inla_spde$param.inla$M1, "G2"=inla_spde$param.inla$M2 )
  if(Version%in%c("spatial_dfa_v13")) TmbData = list("Options_vec"=Options_vec, "n_obs"=Nobs, "n_samples"=Nsamples, "n_obsfactors"=Nobsfactors_input, "n_sites"=Nsites, "n_years"=Nyears, "n_knots"=Nknots, "n_species"=Nspecies, "n_factors"=Nfactors_input, "n_catchcov"=ncol(X_ik), "n_spacecov"=ncol(X_nl), "c_i"=DF[,'catch'], "m_i"=as.numeric(DF[,'TowID'])-1, "p_i"=as.numeric(DF[,'spp'])-1, "s_i"=DF[,'sitenum']-1, "t_i"=match(DF[,'year'],YearSet)-1, "X_ik"=X_ik, "X_nl"=X_nl, "a_n"=a_n, "Rotation_jj"=diag(Nfactors_input), "G0"=inla_spde$param.inla$M0, "G1"=inla_spde$param.inla$M1, "G2"=inla_spde$param.inla$M2 )
  if(Version%in%c("spatial_dfa_v14")) TmbData = list("Options_vec"=Options_vec, "n_obs"=Nobs, "n_samples"=Nsamples, "n_obsfactors"=Nobsfactors_input, "n_sites"=Nsites, "n_years"=Nyears, "n_knots"=Nknots, "n_species"=Nspecies, "n_factors"=Nfactors_input, "n_catchcov"=ncol(X_ik), "n_spacecov"=dim(X_ntl)[3], "c_i"=DF[,'catch'], "m_i"=as.numeric(DF[,'TowID'])-1, "p_i"=as.numeric(DF[,'spp'])-1, "s_i"=DF[,'sitenum']-1, "t_i"=match(DF[,'year'],YearSet)-1, "X_ik"=X_ik, "X_ntl"=X_ntl, "a_n"=a_n, "Rotation_jj"=diag(Nfactors_input), "G0"=inla_spde$param.inla$M0, "G1"=inla_spde$param.inla$M1, "G2"=inla_spde$param.inla$M2 )
  if(Version%in%c("spatial_dfa_v15")) TmbData = list("Options_vec"=Options_vec, "n_obs"=Nobs, "n_samples"=Nsamples, "n_obsfactors"=Nobsfactors_input, "n_sites"=Nsites, "n_years"=Nyears, "n_knots"=Nknots, "n_species"=Nspecies, "n_factors"=Nfactors_input, "n_catchcov"=ncol(X_ik), "n_spacecov"=dim(X_ntl)[3], "c_i"=DF[,'catch'], "predTF_i"=DF[,'PredTF_i'], "m_i"=as.numeric(DF[,'TowID'])-1, "p_i"=as.numeric(DF[,'spp'])-1, "s_i"=DF[,'sitenum']-1, "t_i"=match(DF[,'year'],YearSet)-1, "X_ik"=X_ik, "X_ntl"=X_ntl, "a_n"=a_n, "Rotation_jj"=diag(Nfactors_input), "G0"=inla_spde$param.inla$M0, "G1"=inla_spde$param.inla$M1, "G2"=inla_spde$param.inla$M2 )
  if(Version%in%c("spatial_dfa_v16")) TmbData = list("Options_vec"=Options_vec, "n_obs"=Nobs, "n_samples"=Nsamples, "n_obsfactors"=Nobsfactors_input, "n_sites"=Nsites, "n_years"=Nyears, "n_knots"=Nknots, "n_species"=Nspecies, "n_factors"=Nfactors_input, "n_catchcov"=ncol(X_ik), "n_spacecov"=dim(X_ntl)[3], "c_i"=DF[,'catch'], "predTF_i"=DF[,'PredTF_i'], "m_i"=as.numeric(DF[,'TowID'])-1, "p_i"=as.numeric(DF[,'spp'])-1, "s_i"=DF[,'sitenum']-1, "t_i"=match(DF[,'year'],YearSet)-1, "X_ik"=X_ik, "X_ntl"=X_ntl, "a_n"=a_n, "Rotation_jj"=diag(Nfactors_input), "spde"=NULL, "spde_aniso"=NULL )
  if(Version%in%c("spatial_dfa_v17")) TmbData = list("Options_vec"=Options_vec, "n_obs"=Nobs, "n_samples"=Nsamples, "n_obsfactors"=Nobsfactors_input, "n_sites"=Nsites, "n_years"=Nyears, "n_knots"=Nknots, "n_species"=Nspecies, "n_factors"=Nfactors_input, "n_catchcov"=ncol(X_ik), "n_spacecov"=dim(X_ntl)[3], "c_i"=DF[,'catch'], "predTF_i"=DF[,'PredTF_i'], "m_i"=as.numeric(DF[,'TowID'])-1, "p_i"=as.numeric(DF[,'spp'])-1, "s_i"=DF[,'sitenum']-1, "t_i"=match(DF[,'year'],YearSet)-1, "X_ik"=X_ik, "X_ntl"=X_ntl, "a_n"=a_n, "Rotation_jj"=diag(Nfactors_input), "CorrGroup_pp"=CorrGroup_pp, "spde"=NULL, "spde_aniso"=NULL )
  if(Version%in%c("spatial_dfa_v18")) TmbData = list("Options_vec"=Options_vec, "n_obs"=Nobs, "n_samples"=Nsamples, "n_obsfactors"=Nobsfactors_input, "n_sites"=Nsites, "n_years"=Nyears, "n_knots"=Nknots, "n_species"=Nspecies, "n_factors"=Nfactors_input, "n_catchcov"=ncol(X_ik), "n_spacecov"=dim(X_ntl)[3], "c_i"=DF[,'catch'], "predTF_i"=DF[,'PredTF_i'], "m_i"=as.numeric(DF[,'TowID'])-1, "p_i"=as.numeric(DF[,'spp'])-1, "s_i"=DF[,'sitenum']-1, "t_i"=match(DF[,'year'],YearSet)-1, "X_ik"=X_ik, "X_ntl"=X_ntl, "a_n"=a_n, "Rotation_jj"=diag(Nfactors_input), "CorrGroup_pp"=CorrGroup_pp, "spde"=NULL, "spde_aniso"=NULL, "M0"=grid_list$M0, "M1"=grid_list$M1, "M2"=grid_list$M2 )

  # Add aniso inputs
  if( "spde" %in% names(TmbData)){
    if( !("SpatialDeltaGLMM"%in%installed.packages()) || !("mesh"%in%names(formals(SpatialDeltaGLMM::Calc_Anisotropic_Mesh))) ){
      #if("package:SpatialDeltaGLMM"%in%search()) detach("package:SpatialDeltaGLMM")
      message( "Installing package SpatialDeltaGLMM from github" )
      devtools::install_github("nwfsc-assess/geostatistical_delta-GLMM")
    }
    MeshList = SpatialDeltaGLMM::Calc_Anisotropic_Mesh(loc_x=NA, mesh=inla_mesh)
    TmbData[["spde"]] = INLA::inla.spde2.matern(MeshList$mesh)$param.inla[c("M0","M1","M2")]
    TmbData[["spde_aniso"]] = list("n_s"=MeshList$spde$n.spde, "n_tri"=nrow(MeshList$mesh$graph$tv), "Tri_Area"=MeshList$Tri_Area, "E0"=MeshList$E0, "E1"=MeshList$E1, "E2"=MeshList$E2, "TV"=MeshList$TV-1, "G0"=MeshList$spde$param.inla$M0, "G0_inv"=inla.as.dgTMatrix(solve(MeshList$spde$param.inla$M0)) )
  }

  # Parameters
  if(Version=="spatial_dfa_v2") TmbParams = list("logkappa_j"=rep(log(100),TmbData$n_factors), "alpha_j"=rep(0,TmbData$n_factors), "phi_j"=rep(0,TmbData$n_factors), "loglambda_j"=rep(log(1),TmbData$n_factors), "rho_j"=rep(0.2,TmbData$n_factors), "Psi_val"=rnorm(TmbData$n_factors*TmbData$n_species-TmbData$n_factors*(TmbData$n_factors-1)/2), "gamma_k"=rep(1,TmbData$n_cov), "log_sigma_p"=rep(log(1),TmbData$n_species), "Epsilon_input"=array(0,dim=c(TmbData$n_knots,TmbData$n_factors,TmbData$n_years)), "Omega_input"=array(0,dim=c(TmbData$n_knots,TmbData$n_factors)), "delta_i"=rep(0,TmbData$n_obs) )
  if(Version=="spatial_dfa_v3") TmbParams = list("logkappa_j"=rep(log(100),TmbData$n_factors), "alpha_j"=rep(0,TmbData$n_factors), "phi_j"=rep(0,TmbData$n_factors), "loglambda_j"=rep(log(1),TmbData$n_factors), "rho_j"=rep(0.2,TmbData$n_factors), "L_val"=rnorm(TmbData$n_factors*TmbData$n_species-TmbData$n_factors*(TmbData$n_factors-1)/2), "gamma_k"=rep(1,TmbData$n_cov), "log_sigma_p"=rep(log(1),TmbData$n_species), "Epsilon_input"=array(0,dim=c(TmbData$n_knots,TmbData$n_factors,TmbData$n_years)), "Omega_input"=array(0,dim=c(TmbData$n_knots,TmbData$n_factors)), "delta_i"=rep(0,TmbData$n_obs) )
  if(Version%in%c("spatial_dfa_v7","spatial_dfa_v6","spatial_dfa_v5","spatial_dfa_v4")) TmbParams = list("logkappa_j"=rep(log(100),TmbData$n_factors), "alpha_j"=rep(0,TmbData$n_factors), "phi_j"=rep(0,TmbData$n_factors), "loglambda_j"=rep(log(1),TmbData$n_factors), "rho_j"=rep(0.2,TmbData$n_factors), "L_val"=rnorm(TmbData$n_factors*TmbData$n_species-TmbData$n_factors*(TmbData$n_factors-1)/2), "gamma_k"=rep(1,TmbData$n_cov), "log_sigma_p"=rep(log(1),TmbData$n_species), "zinfl_pz"=matrix(0,nrow=TmbData$n_species,ncol=2), "Epsilon_input"=array(0,dim=c(TmbData$n_knots,TmbData$n_factors,TmbData$n_years)), "Omega_input"=array(0,dim=c(TmbData$n_knots,TmbData$n_factors)), "delta_i"=rep(0,TmbData$n_obs) )
  if(Version%in%c("spatial_dfa_v10","spatial_dfa_v9","spatial_dfa_v8b","spatial_dfa_v8")) TmbParams = list("logkappa_jz"=array(log(100),dim=c(TmbData$n_factors,2)), "alpha_j"=rep(0,TmbData$n_factors), "phi_j"=rep(0,TmbData$n_factors), "loglambda_j"=rep(log(1),TmbData$n_factors), "rho_j"=rep(0.2,TmbData$n_factors), "L_val"=rnorm(TmbData$n_factors*TmbData$n_species-TmbData$n_factors*(TmbData$n_factors-1)/2), "gamma_k"=rep(1,TmbData$n_cov), "log_sigma_p"=rep(log(1),TmbData$n_species), "zinfl_pz"=matrix(0,nrow=TmbData$n_species,ncol=2), "Epsilon_input"=array(0,dim=c(TmbData$n_knots,TmbData$n_factors,TmbData$n_years)), "Omega_input"=array(0,dim=c(TmbData$n_knots,TmbData$n_factors)), "delta_i"=rep(0,TmbData$n_obs) )
  if(Version%in%c("spatial_dfa_v12","spatial_dfa_v11")) TmbParams = list("logkappa_jz"=array(log(100),dim=c(TmbData$n_factors,2)), "alpha_j"=rep(0,TmbData$n_factors), "phi_j"=rep(0,TmbData$n_factors), "loglambda_j"=rep(log(1),TmbData$n_factors), "rho_j"=rep(0.2,TmbData$n_factors), "L_val"=rnorm(TmbData$n_factors*TmbData$n_species-TmbData$n_factors*(TmbData$n_factors-1)/2), "L2_val"=rnorm(TmbData$n_obsfactors*TmbData$n_species-TmbData$n_obsfactors*(TmbData$n_obsfactors-1)/2), "gamma_k"=rep(1,TmbData$n_cov), "log_sigma_p"=rep(log(1),TmbData$n_species), "zinfl_pz"=matrix(0,nrow=TmbData$n_species,ncol=2), "Epsilon_input"=array(0,dim=c(TmbData$n_knots,TmbData$n_factors,TmbData$n_years)), "Omega_input"=array(0,dim=c(TmbData$n_knots,TmbData$n_factors)), "delta_i"=rep(0,TmbData$n_obs), "eta_mb"=array(0,dim=c(TmbData$n_samples,TmbData$n_obsfactors)) )
  if(Version%in%c("spatial_dfa_v13")) TmbParams = list("logkappa_jz"=array(log(100),dim=c(TmbData$n_factors,2)), "alpha_j"=rep(0,TmbData$n_factors), "phi_j"=rep(0,TmbData$n_factors), "loglambda_j"=rep(log(1),TmbData$n_factors), "rho_j"=rep(0.2,TmbData$n_factors), "L_val"=rnorm(TmbData$n_factors*TmbData$n_species-TmbData$n_factors*(TmbData$n_factors-1)/2), "L2_val"=rnorm(TmbData$n_obsfactors*TmbData$n_species-TmbData$n_obsfactors*(TmbData$n_obsfactors-1)/2), "gamma_k"=rep(1,TmbData$n_catchcov), "gamma_lp"=matrix(1,nrow=TmbData$n_spacecov,ncol=TmbData$n_species), "log_sigma_p"=rep(log(1),TmbData$n_species), "zinfl_pz"=matrix(0,nrow=TmbData$n_species,ncol=2), "Epsilon_input"=array(0,dim=c(TmbData$n_knots,TmbData$n_factors,TmbData$n_years)), "Omega_input"=array(0,dim=c(TmbData$n_knots,TmbData$n_factors)), "delta_i"=rep(0,TmbData$n_obs), "eta_mb"=array(0,dim=c(TmbData$n_samples,TmbData$n_obsfactors)) )
  if(Version%in%c("spatial_dfa_v15","spatial_dfa_v14")) TmbParams = list("logkappa_jz"=array(log(100),dim=c(TmbData$n_factors,2)), "alpha_j"=rep(0,TmbData$n_factors), "phi_j"=rep(0,TmbData$n_factors), "loglambda_j"=rep(log(1),TmbData$n_factors), "rho_j"=rep(0.2,TmbData$n_factors), "L_val"=rnorm(TmbData$n_factors*TmbData$n_species-TmbData$n_factors*(TmbData$n_factors-1)/2), "L2_val"=rnorm(TmbData$n_obsfactors*TmbData$n_species-TmbData$n_obsfactors*(TmbData$n_obsfactors-1)/2), "gamma_k"=rep(1,TmbData$n_catchcov), "gamma_ptl"=array(1,dim=unlist(TmbData[c("n_species","n_years","n_spacecov")])), "log_sigma_p"=rep(log(1),TmbData$n_species), "zinfl_pz"=matrix(0,nrow=TmbData$n_species,ncol=2), "Epsilon_input"=array(0,dim=c(TmbData$n_knots,TmbData$n_factors,TmbData$n_years)), "Omega_input"=array(0,dim=c(TmbData$n_knots,TmbData$n_factors)), "delta_i"=rep(0,TmbData$n_obs), "eta_mb"=array(0,dim=c(TmbData$n_samples,TmbData$n_obsfactors)) )
  if(Version%in%c("spatial_dfa_v18","spatial_dfa_v17","spatial_dfa_v16")) TmbParams = list("ln_H_input"=c(0,0), "logkappa_jz"=array(-1,dim=c(TmbData$n_factors,2)), "alpha_j"=rep(0,TmbData$n_factors), "phi_j"=rep(0,TmbData$n_factors), "loglambda_j"=rep(log(1),TmbData$n_factors), "rho_j"=rep(0.2,TmbData$n_factors), "L_val"=rnorm(TmbData$n_factors*TmbData$n_species-TmbData$n_factors*(TmbData$n_factors-1)/2), "L2_val"=rnorm(TmbData$n_obsfactors*TmbData$n_species-TmbData$n_obsfactors*(TmbData$n_obsfactors-1)/2), "gamma_k"=rep(1,TmbData$n_catchcov), "gamma_ptl"=array(1,dim=unlist(TmbData[c("n_species","n_years","n_spacecov")])), "log_sigma_p"=rep(log(1),TmbData$n_species), "zinfl_pz"=matrix(0,nrow=TmbData$n_species,ncol=2), "Epsilon_input"=array(0,dim=c(TmbData$n_knots,TmbData$n_factors,TmbData$n_years)), "Omega_input"=array(0,dim=c(TmbData$n_knots,TmbData$n_factors)), "delta_i"=rep(0,TmbData$n_obs), "eta_mb"=array(0,dim=c(TmbData$n_samples,TmbData$n_obsfactors)) )
                                                                                                                                                                   #   ,
  # Random
  Random = c( "Omega_input", "Epsilon_input", "delta_i" )
  if(Use_REML==TRUE) Random = c(Random, "gamma_k")  #  , "log_zinfl"
  if("eta_mb"%in%names(TmbParams)) Random = c(Random, "eta_mb")
  if("gamma_lp"%in%names(TmbParams) && Use_REML==TRUE) Random = c(Random, "gamma_lp")  #  , "log_zinfl"
  if("gamma_ptl"%in%names(TmbParams) && Use_REML==TRUE) Random = c(Random, "gamma_ptl")  #  , "log_zinfl"
  #Random = NULL

  # Fixed values
  Map = list()
  Map[["alpha_j"]] = factor( rep(NA,length(TmbParams[["alpha_j"]])) )
  #Map[["phi_j"]] = factor( rep(NA,length(TmbParams[["phi_j"]])) )
  # Parameters shared between dynamic factors
  Map[["rho_j"]] = factor( rep(1,length(TmbParams[["rho_j"]])) )
  Map[["loglambda_j"]] = factor( rep(1,length(TmbParams[["loglambda_j"]])) )
  if(Version%in%c("spatial_dfa_v7","spatial_dfa_v6","spatial_dfa_v5","spatial_dfa_v4")){
    if( Kappa_Type=="Constant" ) Map[["logkappa_j"]] = factor( rep(1,length(TmbParams[["logkappa_j"]])) )
    if( Kappa_Type=="Omega_vs_Epsilon" ) stop("Not implemented")
  }
  if(Version%in%c("spatial_dfa_v18","spatial_dfa_v17","spatial_dfa_v16","spatial_dfa_v15","spatial_dfa_v14","spatial_dfa_v13","spatial_dfa_v12","spatial_dfa_v11","spatial_dfa_v10","spatial_dfa_v9","spatial_dfa_v8b","spatial_dfa_v8")){
    if( Kappa_Type=="Constant" ) Map[["logkappa_jz"]] = factor( array(1,dim=dim(TmbParams[["logkappa_jz"]])) )
    if( Kappa_Type=="Omega_vs_Epsilon" ) Map[["logkappa_jz"]] = factor( outer(rep(1,TmbData$n_factors),c(1,2)) )
  }

  ##### Run-specific fixed values
  # Turn off anisotropy
  if( "ln_H_input"%in%names(TmbParams) & Aniso==FALSE ){
    Map[["ln_H_input"]] = factor( c(NA,NA) )
  }
  # Turn off zero-inflation parameters
  if( !is.na(Options_vec["EncounterFunction"]) && Options_vec["EncounterFunction"]==2 ){
    Map[["zinfl_pz"]] = factor( cbind(1:TmbData$n_species, rep(NA,TmbData$n_species)) )
  }
  # ObsModel specification
  if( Options_vec['ObsModel'] %in% c(0,4) ){
    Map[["zinfl_pz"]] = factor( rep(NA,prod(dim(TmbParams[["zinfl_pz"]]))) )
  }
  if( Options_vec['ObsModel'] %in% c(1,2,3,4) ){
    # Shrink size for speed-up during compile
    TmbParams[["delta_i"]] = 0
    Map[["delta_i"]] = factor( rep(NA,length(TmbParams[["delta_i"]])) )
  }
  # Turn off Phi
  if( Include_Phi==FALSE | Nfactors==0 ){
    Map[["phi_j"]] = factor( rep(NA,length(TmbParams[["phi_j"]])) )
    TmbParams[["phi_j"]] = rep(0,length(TmbParams[["phi_j"]]))
  }
  # Turn off Rho
  if( Include_Rho==FALSE | Nfactors==0 ){
    Map[["rho_j"]] = factor( rep(NA,length(TmbParams[["rho_j"]])) )
    TmbParams[["rho_j"]] = rep(0,length(TmbParams[["rho_j"]]))
  }
  # Turn off spatial variation (Omega)
  if( Options_vec['Include_Omega']==FALSE | Nfactors==0 ){
    Map[["Omega_input"]] = factor( array(NA,dim=dim(TmbParams[["Omega_input"]])) )
    TmbParams[["Omega_input"]] = array(0,dim=dim(TmbParams[["Omega_input"]]))
    Map[["loglambda_j"]] = factor( rep(NA,length(TmbParams[["loglambda_j"]])) )
    TmbParams[["loglambda_j"]] = rep(25,length(TmbParams[["loglambda_j"]]))
  }
  # Turn off spatiotemporal variation (Epsilon)
  if( Options_vec['Include_Epsilon']==FALSE | Nfactors==0 ){
    Map[["Epsilon_input"]] = factor( array(NA,dim=dim(TmbParams[["Epsilon_input"]])) )
    TmbParams[["Epsilon_input"]] = array(0,dim=dim(TmbParams[["Epsilon_input"]]))
    Map[["loglambda_j"]] = factor( rep(NA,length(TmbParams[["loglambda_j"]])) )
    TmbParams[["loglambda_j"]] = rep(25,length(TmbParams[["loglambda_j"]]))
  }
  # Turn off both spatial and spatiotemporal variation
  if( (Options_vec['Include_Omega']==FALSE & Options_vec['Include_Epsilon']==FALSE) | Nfactors==0 ){
    TmbData$Options_vec[c('Include_Omega','Include_Epsilon')] = 0
    Map[["logkappa_jz"]] = factor( array(NA,dim=dim(TmbParams[["logkappa_jz"]])) )
    Map[["rho_j"]] = factor( rep(NA,length(TmbParams[["rho_j"]])) )
    Map[["L_val"]] = factor( rep(NA,length(TmbParams[["L_val"]])) )
    if( "ln_H_input"%in%names(TmbParams) ) Map[["ln_H_input"]] = factor( c(NA,NA) )
  }
  # 
  if( Options_vec["Correlated_Overdispersion"]==0 ){
    TmbParams[["L2_val"]][] = 0
    Map[["L2_val"]] = factor( rep(NA,length(TmbParams[["L2_val"]])) )
    # Shrink size for speed-up during compile
    TmbParams[["eta_mb"]] = array(0,dim=c(1,ncol(TmbParams$eta_mb)))
    Map[["eta_mb"]] = factor( array(NA,dim(TmbParams[["eta_mb"]])) )
  }
  # Default setting for spatio-temporal covariates -- constant across time
  if( "gamma_ptl"%in%names(TmbParams) ){
    Map[["gamma_ptl"]] = factor( aperm(outer(matrix(1:(TmbData$n_species*TmbData$n_spacecov),ncol=TmbData$n_spacecov,nrow=TmbData$n_species),rep(1,TmbData$n_years)),c(1,3,2)) )
  }

  # Turn off catchability covariates
  if( "gamma_k"%in%names(TmbParams) && all(X_ik==0) ){
    Map[["gamma_k"]] = factor( rep(NA,length(TmbParams[["gamma_k"]])) )
    TmbParams[["gamma_k"]][] = 0
  }
  # Turn off spatial covariates
  if( "gamma_lp"%in%names(TmbParams) && all(X_nl==0) ){
    Map[["gamma_lp"]] = factor( array(NA,dim=dim(TmbParams[["gamma_lp"]])) )
    TmbParams[["gamma_lp"]][] = 0
  }
  # Turn off spatio-temporal covariates
  if( "gamma_ptl"%in%names(TmbParams) && all(X_ntl==0) ){
    Map[["gamma_ptl"]] = factor( array(NA,dim=dim(TmbParams[["gamma_ptl"]])) )
    TmbParams[["gamma_ptl"]][] = 0
  }

  # Make independent inputs
  if( IndependentTF[1]==TRUE ){
    if(Nfactors!=Nspecies) stop("If independent, Nfactors must equal Nspecies")
    TmbParams[["L_val"]] = diag( rep(1,Nspecies) )[lower.tri(diag(1,Nspecies),diag=TRUE)]
    Map[["L_val"]] = diag( 1:Nspecies )[lower.tri(diag(1,Nspecies),diag=TRUE)]
    Map[["L_val"]] = factor( ifelse(Map[["L_val"]]==0,NA,Map[["L_val"]]) )
  }
  if( IndependentTF[2]==TRUE ){
    if(Nobsfactors!=Nspecies) stop("If independent, Nfactors must equal Nspecies")
    TmbParams[["L2_val"]] = diag( rep(1,Nspecies) )[lower.tri(diag(1,Nspecies),diag=TRUE)]
    Map[["L2_val"]] = diag( 1:Nspecies )[lower.tri(diag(1,Nspecies),diag=TRUE)]
    Map[["L2_val"]] = factor( ifelse(Map[["L_val"]]==0,NA,Map[["L_val"]]) )
  }

  # Check for bugs
  if( CheckForBugs==TRUE ){
    if( any(sapply(TmbParams[c("alpha_j","phi_j","loglambda_j","rho_j")],length)<TmbData$n_factors) ) stop("Problem with parameter-vectors subscripted j")
    if( is.null(EncounterFunction) | is.null(ObsModel)) stop("Problem with NULL inputs")
    if( Options_vec["Correlated_Overdispersion"]==1 ) if( max(TmbData$m_i)>TmbData$n_samples | min(TmbData$m_i)<0 ) stop("Problem with m_i")
    if( max(TmbData$t_i)>TmbData$n_years | min(TmbData$t_i)<0 ) stop("Problem with t_i")
  }

  # Check mapped stuff
  Remove_Random = NULL
  for(i in 1:length(Random) ){
    if( Random[i]%in%names(Map) && all(is.na(Map[[Random[i]]])) ){
      Remove_Random = c(Remove_Random, Random[i])
    }
  }
  Random = setdiff(Random, Remove_Random)
  if(length(Random)==0) Random = NULL

  # Return
  Return = list("TmbData"=TmbData, "TmbParams"=TmbParams, "Random"=Random, "Map"=Map, "mesh"=inla_mesh, "Remove_Random"=Remove_Random)
  return( Return )
}

