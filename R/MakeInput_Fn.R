MakeInput_Fn = function( Version, Nfactors, Nobsfactors=0, DF, inla_spde, Kappa_Type="Constant", ObsModel=NULL, Include_Omega=TRUE, Include_Epsilon=TRUE, EncounterFunction=2, Correlated_Overdispersion=FALSE, Include_Phi=TRUE, Include_Rho=TRUE, Use_REML=FALSE, X_ik=NULL, X_nl=NULL, X_ntl=NULL, a_n=NULL, YearSet=NULL, CheckForBugs=TRUE){

  # Options_vec
  Options_vec=c( "ObsModel"=ObsModel, "Include_Omega"=Include_Omega, "Include_Epsilon"=Include_Epsilon, "EncounterFunction"=EncounterFunction, "Correlated_Overdispersion"=ifelse(Nobsfactors==0,0,1))

  # Infer default values for inputs
  if( is.null(YearSet) ) YearSet = min(DF[,'year']):max(DF[,'year'])
  if( is.null(ObsModel) ){
    ObsModel = ifelse( all(is.integer(ObsModel[,'catch'])), 0, 1 )
  }
  if( is.null(a_n) ) a_n = rep(0,mesh$n)

  # Data size
  Nyears = length(YearSet)
  Nsites = length(unique(DF[,'sitenum']))
  Nspecies = length(levels(DF[,'spp']))
  Nknots = mesh$n
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
  
  # Default catchability design matrix
  if( is.null(X_ik) ){
    X_ik = matrix(0, nrow=Nobs, ncol=1)
  }else{
    if( nrow(X_ik)!=Nobs ) stop("Check catchability design matrix input: X_ik")
  }

  # Default spatial design matrix
  if( is.null(X_nl) ){
    X_nl = matrix(1, nrow=Nknots, ncol=1)
  }else{
    if( !(Version%in%c("spatial_dfa_v13")) ) stop("X_nl is only used in version spatial_dfa_v13")
    if( nrow(X_nl)!=Nknots ) stop("Check spatial matrix input: X_sp")    
  }

  # Default spatial design matrix
  if( is.null(X_ntl) ){
    X_ntl = array(1, dim=c(Nknots,Nyears,1))
  }else{
    if( !(Version%in%c("spatial_dfa_v14")) ) stop("X_ntl is only used in version spatial_dfa_v14 and higher")
    if( dim(X_ntl)[1] != Nknots ) stop("Check spatial matrix input: X_ntl")    
  }

  # Data
  if(Version=="spatial_dfa_v3") TmbData = list("n_obs"=Nobs, "n_sites"=Nsites, "n_years"=Nyears, "n_knots"=Nknots, "n_species"=Nspecies, "n_factors"=Nfactors_input, "n_cov"=ncol(X_ik), "c_i"=DF[,'catch'], "p_i"=as.numeric(DF[,'spp'])-1, "s_i"=DF[,'sitenum']-1, "t_i"=match(DF[,'year'],YearSet)-1, "X_ik"=X_ik, "G0"=inla_spde$param.inla$M0, "G1"=inla_spde$param.inla$M1, "G2"=inla_spde$param.inla$M2 )
  if(Version=="spatial_dfa_v4") TmbData = list("Options_vec"=Options_vec, "n_obs"=Nobs, "n_sites"=Nsites, "n_years"=Nyears, "n_knots"=Nknots, "n_species"=Nspecies, "n_factors"=Nfactors_input, "n_cov"=ncol(X_ik), "c_i"=DF[,'catch'], "p_i"=as.numeric(DF[,'spp'])-1, "s_i"=DF[,'sitenum']-1, "t_i"=match(DF[,'year'],YearSet)-1, "X_ik"=X_ik, "G0"=inla_spde$param.inla$M0, "G1"=inla_spde$param.inla$M1, "G2"=inla_spde$param.inla$M2 )
  if(Version%in%c("spatial_dfa_v10","spatial_dfa_v9","spatial_dfa_v8b","spatial_dfa_v8","spatial_dfa_v7","spatial_dfa_v6","spatial_dfa_v5")) TmbData = list("Options_vec"=Options_vec, "n_obs"=Nobs, "n_sites"=Nsites, "n_years"=Nyears, "n_knots"=Nknots, "n_species"=Nspecies, "n_factors"=Nfactors_input, "n_cov"=ncol(X_ik), "c_i"=DF[,'catch'], "p_i"=as.numeric(DF[,'spp'])-1, "s_i"=DF[,'sitenum']-1, "t_i"=match(DF[,'year'],YearSet)-1, "X_ik"=X_ik, "Rotation_jj"=diag(Nfactors_input), "G0"=inla_spde$param.inla$M0, "G1"=inla_spde$param.inla$M1, "G2"=inla_spde$param.inla$M2 )
  if(Version%in%c("spatial_dfa_v11")) TmbData = list("Options_vec"=Options_vec, "n_obs"=Nobs, "n_samples"=Nsamples, "n_obsfactors"=Nobsfactors_input, "n_sites"=Nsites, "n_years"=Nyears, "n_knots"=Nknots, "n_species"=Nspecies, "n_factors"=Nfactors_input, "n_cov"=ncol(X_ik), "c_i"=DF[,'catch'], "m_i"=as.numeric(DF[,'TowID'])-1, "p_i"=as.numeric(DF[,'spp'])-1, "s_i"=DF[,'sitenum']-1, "t_i"=match(DF[,'year'],YearSet)-1, "X_ik"=X_ik, "Rotation_jj"=diag(Nfactors_input), "G0"=inla_spde$param.inla$M0, "G1"=inla_spde$param.inla$M1, "G2"=inla_spde$param.inla$M2 )
  if(Version%in%c("spatial_dfa_v12")) TmbData = list("Options_vec"=Options_vec, "n_obs"=Nobs, "n_samples"=Nsamples, "n_obsfactors"=Nobsfactors_input, "n_sites"=Nsites, "n_years"=Nyears, "n_knots"=Nknots, "n_species"=Nspecies, "n_factors"=Nfactors_input, "n_cov"=ncol(X_ik), "c_i"=DF[,'catch'], "m_i"=as.numeric(DF[,'TowID'])-1, "p_i"=as.numeric(DF[,'spp'])-1, "s_i"=DF[,'sitenum']-1, "t_i"=match(DF[,'year'],YearSet)-1, "X_ik"=X_ik, "a_n"=a_n, "Rotation_jj"=diag(Nfactors_input), "G0"=inla_spde$param.inla$M0, "G1"=inla_spde$param.inla$M1, "G2"=inla_spde$param.inla$M2 )
  if(Version%in%c("spatial_dfa_v13")) TmbData = list("Options_vec"=Options_vec, "n_obs"=Nobs, "n_samples"=Nsamples, "n_obsfactors"=Nobsfactors_input, "n_sites"=Nsites, "n_years"=Nyears, "n_knots"=Nknots, "n_species"=Nspecies, "n_factors"=Nfactors_input, "n_catchcov"=ncol(X_ik), "n_spacecov"=ncol(X_nl), "c_i"=DF[,'catch'], "m_i"=as.numeric(DF[,'TowID'])-1, "p_i"=as.numeric(DF[,'spp'])-1, "s_i"=DF[,'sitenum']-1, "t_i"=match(DF[,'year'],YearSet)-1, "X_ik"=X_ik, "X_nl"=X_nl, "a_n"=a_n, "Rotation_jj"=diag(Nfactors_input), "G0"=inla_spde$param.inla$M0, "G1"=inla_spde$param.inla$M1, "G2"=inla_spde$param.inla$M2 )
  if(Version%in%c("spatial_dfa_v14")) TmbData = list("Options_vec"=Options_vec, "n_obs"=Nobs, "n_samples"=Nsamples, "n_obsfactors"=Nobsfactors_input, "n_sites"=Nsites, "n_years"=Nyears, "n_knots"=Nknots, "n_species"=Nspecies, "n_factors"=Nfactors_input, "n_catchcov"=ncol(X_ik), "n_spacecov"=ncol(X_nl), "c_i"=DF[,'catch'], "m_i"=as.numeric(DF[,'TowID'])-1, "p_i"=as.numeric(DF[,'spp'])-1, "s_i"=DF[,'sitenum']-1, "t_i"=match(DF[,'year'],YearSet)-1, "X_ik"=X_ik, "X_ntl"=X_ntl, "a_n"=a_n, "Rotation_jj"=diag(Nfactors_input), "G0"=inla_spde$param.inla$M0, "G1"=inla_spde$param.inla$M1, "G2"=inla_spde$param.inla$M2 )
  #TmbData$c_i = c(1, rep(NA,TmbData$n_obs-1) )

  # Parameters
  if(Version=="spatial_dfa_v2") TmbParams = list("logkappa_j"=rep(log(100),TmbData$n_factors), "alpha_j"=rep(0,TmbData$n_factors), "phi_j"=rep(0,TmbData$n_factors), "loglambda_j"=rep(log(1),TmbData$n_factors), "rho_j"=rep(0.2,TmbData$n_factors), "Psi_val"=rnorm(TmbData$n_factors*TmbData$n_species-TmbData$n_factors*(TmbData$n_factors-1)/2), "gamma_k"=rep(1,TmbData$n_cov), "log_sigma_p"=rep(log(1),TmbData$n_species), "Epsilon_input"=array(0,dim=c(TmbData$n_knots,TmbData$n_factors,TmbData$n_years)), "Omega_input"=array(0,dim=c(TmbData$n_knots,TmbData$n_factors)), "delta_i"=rep(0,TmbData$n_obs) )
  if(Version=="spatial_dfa_v3") TmbParams = list("logkappa_j"=rep(log(100),TmbData$n_factors), "alpha_j"=rep(0,TmbData$n_factors), "phi_j"=rep(0,TmbData$n_factors), "loglambda_j"=rep(log(1),TmbData$n_factors), "rho_j"=rep(0.2,TmbData$n_factors), "L_val"=rnorm(TmbData$n_factors*TmbData$n_species-TmbData$n_factors*(TmbData$n_factors-1)/2), "gamma_k"=rep(1,TmbData$n_cov), "log_sigma_p"=rep(log(1),TmbData$n_species), "Epsilon_input"=array(0,dim=c(TmbData$n_knots,TmbData$n_factors,TmbData$n_years)), "Omega_input"=array(0,dim=c(TmbData$n_knots,TmbData$n_factors)), "delta_i"=rep(0,TmbData$n_obs) )
  if(Version%in%c("spatial_dfa_v7","spatial_dfa_v6","spatial_dfa_v5","spatial_dfa_v4")) TmbParams = list("logkappa_j"=rep(log(100),TmbData$n_factors), "alpha_j"=rep(0,TmbData$n_factors), "phi_j"=rep(0,TmbData$n_factors), "loglambda_j"=rep(log(1),TmbData$n_factors), "rho_j"=rep(0.2,TmbData$n_factors), "L_val"=rnorm(TmbData$n_factors*TmbData$n_species-TmbData$n_factors*(TmbData$n_factors-1)/2), "gamma_k"=rep(1,TmbData$n_cov), "log_sigma_p"=rep(log(1),TmbData$n_species), "zinfl_pz"=matrix(0,nrow=TmbData$n_species,ncol=2), "Epsilon_input"=array(0,dim=c(TmbData$n_knots,TmbData$n_factors,TmbData$n_years)), "Omega_input"=array(0,dim=c(TmbData$n_knots,TmbData$n_factors)), "delta_i"=rep(0,TmbData$n_obs) )
  if(Version%in%c("spatial_dfa_v10","spatial_dfa_v9","spatial_dfa_v8b","spatial_dfa_v8")) TmbParams = list("logkappa_jz"=array(log(100),dim=c(TmbData$n_factors,2)), "alpha_j"=rep(0,TmbData$n_factors), "phi_j"=rep(0,TmbData$n_factors), "loglambda_j"=rep(log(1),TmbData$n_factors), "rho_j"=rep(0.2,TmbData$n_factors), "L_val"=rnorm(TmbData$n_factors*TmbData$n_species-TmbData$n_factors*(TmbData$n_factors-1)/2), "gamma_k"=rep(1,TmbData$n_cov), "log_sigma_p"=rep(log(1),TmbData$n_species), "zinfl_pz"=matrix(0,nrow=TmbData$n_species,ncol=2), "Epsilon_input"=array(0,dim=c(TmbData$n_knots,TmbData$n_factors,TmbData$n_years)), "Omega_input"=array(0,dim=c(TmbData$n_knots,TmbData$n_factors)), "delta_i"=rep(0,TmbData$n_obs) )
  if(Version%in%c("spatial_dfa_v12","spatial_dfa_v11")) TmbParams = list("logkappa_jz"=array(log(100),dim=c(TmbData$n_factors,2)), "alpha_j"=rep(0,TmbData$n_factors), "phi_j"=rep(0,TmbData$n_factors), "loglambda_j"=rep(log(1),TmbData$n_factors), "rho_j"=rep(0.2,TmbData$n_factors), "L_val"=rnorm(TmbData$n_factors*TmbData$n_species-TmbData$n_factors*(TmbData$n_factors-1)/2), "L2_val"=rnorm(TmbData$n_obsfactors*TmbData$n_species-TmbData$n_obsfactors*(TmbData$n_obsfactors-1)/2), "gamma_k"=rep(1,TmbData$n_cov), "log_sigma_p"=rep(log(1),TmbData$n_species), "zinfl_pz"=matrix(0,nrow=TmbData$n_species,ncol=2), "Epsilon_input"=array(0,dim=c(TmbData$n_knots,TmbData$n_factors,TmbData$n_years)), "Omega_input"=array(0,dim=c(TmbData$n_knots,TmbData$n_factors)), "delta_i"=rep(0,TmbData$n_obs), "eta_mb"=array(0,dim=c(TmbData$n_samples,TmbData$n_obsfactors)) )
  if(Version%in%c("spatial_dfa_v13")) TmbParams = list("logkappa_jz"=array(log(100),dim=c(TmbData$n_factors,2)), "alpha_j"=rep(0,TmbData$n_factors), "phi_j"=rep(0,TmbData$n_factors), "loglambda_j"=rep(log(1),TmbData$n_factors), "rho_j"=rep(0.2,TmbData$n_factors), "L_val"=rnorm(TmbData$n_factors*TmbData$n_species-TmbData$n_factors*(TmbData$n_factors-1)/2), "L2_val"=rnorm(TmbData$n_obsfactors*TmbData$n_species-TmbData$n_obsfactors*(TmbData$n_obsfactors-1)/2), "gamma_k"=rep(1,TmbData$n_catchcov), "gamma_lp"=matrix(1,nrow=TmbData$n_spacecov,ncol=TmbData$n_species), "log_sigma_p"=rep(log(1),TmbData$n_species), "zinfl_pz"=matrix(0,nrow=TmbData$n_species,ncol=2), "Epsilon_input"=array(0,dim=c(TmbData$n_knots,TmbData$n_factors,TmbData$n_years)), "Omega_input"=array(0,dim=c(TmbData$n_knots,TmbData$n_factors)), "delta_i"=rep(0,TmbData$n_obs), "eta_mb"=array(0,dim=c(TmbData$n_samples,TmbData$n_obsfactors)) )
  if(Version%in%c("spatial_dfa_v14")) TmbParams = list("logkappa_jz"=array(log(100),dim=c(TmbData$n_factors,2)), "alpha_j"=rep(0,TmbData$n_factors), "phi_j"=rep(0,TmbData$n_factors), "loglambda_j"=rep(log(1),TmbData$n_factors), "rho_j"=rep(0.2,TmbData$n_factors), "L_val"=rnorm(TmbData$n_factors*TmbData$n_species-TmbData$n_factors*(TmbData$n_factors-1)/2), "L2_val"=rnorm(TmbData$n_obsfactors*TmbData$n_species-TmbData$n_obsfactors*(TmbData$n_obsfactors-1)/2), "gamma_k"=rep(1,TmbData$n_catchcov), "gamma_ptl"=array(1,dim=unlist(TmbData[c("n_species","n_years","n_spacecov")])), "log_sigma_p"=rep(log(1),TmbData$n_species), "zinfl_pz"=matrix(0,nrow=TmbData$n_species,ncol=2), "Epsilon_input"=array(0,dim=c(TmbData$n_knots,TmbData$n_factors,TmbData$n_years)), "Omega_input"=array(0,dim=c(TmbData$n_knots,TmbData$n_factors)), "delta_i"=rep(0,TmbData$n_obs), "eta_mb"=array(0,dim=c(TmbData$n_samples,TmbData$n_obsfactors)) )
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
  if(Version%in%c("spatial_dfa_v14","spatial_dfa_v13","spatial_dfa_v12","spatial_dfa_v11","spatial_dfa_v10","spatial_dfa_v9","spatial_dfa_v8b","spatial_dfa_v8")){
    if( Kappa_Type=="Constant" ) Map[["logkappa_jz"]] = factor( array(1,dim=dim(TmbParams[["logkappa_jz"]])) )
    if( Kappa_Type=="Omega_vs_Epsilon" ) Map[["logkappa_jz"]] = factor( outer(rep(1,TmbData$n_factors),c(1,2)) )
  }

  ##### Run-specific fixed values
  # ObsModel specification
  if( Options_vec['ObsModel']==0){
    Map[["zinfl_pz"]] = factor( rep(NA,prod(dim(TmbParams[["zinfl_pz"]]))) )
  }
  if( Options_vec['ObsModel']%in%c(1,2,3) ){
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
  }
  # Turn off zero-inflation parameters
  if( !is.na(Options_vec["EncounterFunction"]) && Options_vec["EncounterFunction"]==2 ){
    Map[["zinfl_pz"]] = factor( cbind(1:TmbData$n_species, rep(NA,TmbData$n_species)) )
  }
  # 
  if( Options_vec["Correlated_Overdispersion"]==0 ){
    TmbParams[["L2_val"]][] = 0
    Map[["L2_val"]] = factor( rep(NA,length(TmbParams[["L2_val"]])) )
    # Shrink size for speed-up during compile
    TmbParams[["eta_mb"]] = array(0,dim=c(1,ncol(TmbParams$eta_mb)))
    Map[["eta_mb"]] = factor( array(NA,dim(TmbParams[["eta_mb"]])) )
  }
  # Default setting for spatial covariates -- constant across time
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
  # Turn off spatial covariates
  if( "gamma_ptl"%in%names(TmbParams) && all(X_ntl==0) ){
    Map[["gamma_ptl"]] = factor( array(NA,dim=dim(TmbParams[["gamma_ptl"]])) )
    TmbParams[["gamma_ptl"]][] = 0
  }

  # Check for bugs
  if( CheckForBugs==TRUE ){
    if( any(sapply(TmbParams[c("alpha_j","phi_j","loglambda_j","rho_j")],length)<TmbData$n_factors) ) stop("Problem with parameter-vectors subscripted j")
    if( is.null(EncounterFunction) | is.null(ObsModel)) stop("Problem with NULL inputs")
    if( max(TmbData$m_i)>TmbData$n_samples | min(TmbData$m_i)<0 ) stop("Problem with m_i")
    if( max(TmbData$t_i)>TmbData$n_years | min(TmbData$t_i)<0 ) stop("Problem with t_i")
  }

  # Return
  Return = list("TmbData"=TmbData, "TmbParams"=TmbParams, "Random"=Random, "Map"=Map)
  return( Return )
}

