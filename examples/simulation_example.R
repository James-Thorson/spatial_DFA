
# Install package
devtools::install_github("james-thorson/Spatial_DFA", auth_token="918291a743d6fe53aecec3bfc4a27d2b850c688d")

# File structure
RootFile = paste0(getwd(),"/")
TmbFile = system.file("executables", package="SpatialDFA")

# Settings
MeshType = c("Minimal", "Recommended")[2]
Version = "spatial_dfa_v4"
  # v3b_failed -- added mu_i
  # v4 -- added zero-inflated version
DataSet = c("Aleknagik", "Bering_Sea", "Simulated")[3]
Kmeans_Config = list("nstart"=10, "iter.max"=1e3, "n_x"=25)
Sim_Settings = list("n_species"=5, "n_years"=20, "n_stations"=20, "n_factors"=2, "SpatialScale"=0.25, "SD_extra"=0.05)

# Model-related settings
Nfactors = 2
Use_REML = TRUE
Include_Temp = FALSE

# Derived settings
ObsModel = switch( DataSet, "Aleknagik"=0, "Bering_Sea"=1, "Simulated"=0)
Nreps = switch( DataSet, "Aleknagik"=1, "Bering_Sea"=1, "Simulated"=100)

# Date file
Date = Sys.Date()
  #Date = "2015-03-17"
  DateFile = paste0(RootFile,Date,"_DataSet=",DataSet,"_IncludeTemp=",Include_Temp,"/")
  dir.create( DateFile )

# Libraries
library( INLA )
library( TMB )

# Specific libraries
library( SpatialDFA )  

# Compile TMB model
setwd( TmbFile )
if(FALSE){
  dyn.unload(dynlib(Version))
  file.remove( paste0(Version,c(".o",".dll")) )
}
compile( paste0(Version,".cpp") )

##############
# Simulation
##############

# Simulate data
Sim_List = Sim_Fn(n_species=Sim_Settings[["n_species"]], n_years=Sim_Settings[["n_years"]], n_stations=Sim_Settings[["n_stations"]], n_factors=Sim_Settings[["n_factors"]], SpatialScale=Sim_Settings[["SpatialScale"]], SD_extra=Sim_Settings[["SD_extra"]])

# Unload stuff
DF = Sim_List[["DF"]]

# lat_set
latlong_set = paste( Sim_List[["Loc"]][,'y'], Sim_List[["Loc"]][,'x'], sep="_")
lat_set = as.numeric(sapply( latlong_set, FUN=function(Char){strsplit(Char,"_")[[1]][1]}))
long_set = as.numeric(sapply( latlong_set, FUN=function(Char){strsplit(Char,"_")[[1]][2]}))

## Build SPDE object using INLA
if(MeshType=="Recommended") mesh = inla.mesh.create( cbind(long_set, lat_set), plot.delay=NULL, extend=list(n=8,offset=-0.15), refine=list(min.angle=26) )  # loc_samp  ;  ,max.edge.data=0.08,max.edge.extra=0.2
if(MeshType=="Minimal") mesh = inla.mesh.create( cbind(long_set, lat_set), plot.delay=NULL, extend=list(n=8,offset=-0.15), refine=F )  # loc_samp

## Create the SPDE/GMRF model, (kappa^2-Delta)(tau x) = W:
spde = inla.spde2.matern(mesh,alpha=2)

# Dimensions
YearSet = min(DF[,'year']):max(DF[,'year'])
Nyears = length(YearSet)
Nsites = length(unique(DF[,'sitenum']))
Nspecies = length(levels(DF[,'spp']))
Nknots = mesh$n
Nobs = nrow(DF)

#### Generate inputs for TMB
# Design matrix
X_ik = matrix(0, nrow=Nobs, ncol=Nspecies)
X_ik[ cbind(1:Nobs,as.numeric(DF[,'spp'])) ] = 1
if( Include_Temp==TRUE ){
  Y_ik = X_ik * outer( DF[,'waterTmpC']-mean(DF[,'waterTmpC']), rep(1,Nspecies) )
  X_ik = cbind( X_ik, Y_ik )
}

# Data
if(Version=="spatial_dfa_v3") TmbData = list("n_obs"=Nobs, "n_sites"=Nsites, "n_years"=Nyears, "n_knots"=Nknots, "n_species"=Nspecies, "n_factors"=Nfactors, "n_cov"=ncol(X_ik), "c_i"=DF[,'catch'], "p_i"=as.numeric(DF[,'spp'])-1, "s_i"=DF[,'sitenum']-1, "t_i"=DF[,'year']-min(DF[,'year']), "X_ik"=X_ik, "G0"=spde$param.inla$M0, "G1"=spde$param.inla$M1, "G2"=spde$param.inla$M2 ) 
if(Version=="spatial_dfa_v4") TmbData = list("Options_vec"=c("ObsModel"=ObsModel), "n_obs"=Nobs, "n_sites"=Nsites, "n_years"=Nyears, "n_knots"=Nknots, "n_species"=Nspecies, "n_factors"=Nfactors, "n_cov"=ncol(X_ik), "c_i"=DF[,'catch'], "p_i"=as.numeric(DF[,'spp'])-1, "s_i"=DF[,'sitenum']-1, "t_i"=DF[,'year']-min(DF[,'year']), "X_ik"=X_ik, "G0"=spde$param.inla$M0, "G1"=spde$param.inla$M1, "G2"=spde$param.inla$M2 ) 

# Parameters
if(Version=="spatial_dfa_v2") TmbParams = list("logkappa_j"=rep(log(100),TmbData$n_factors), "alpha_j"=rep(0,TmbData$n_factors), "phi_j"=rep(0,TmbData$n_factors), "loglambda_j"=rep(log(1),TmbData$n_factors), "rho_j"=rep(0.2,TmbData$n_factors), "Psi_val"=rnorm(TmbData$n_factors*TmbData$n_species-TmbData$n_factors*(TmbData$n_factors-1)/2), "gamma_k"=rep(1,TmbData$n_cov), "log_sigma_p"=rep(log(1),TmbData$n_species), "Epsilon_input"=array(0,dim=c(TmbData$n_knots,TmbData$n_factors,TmbData$n_years)), "Omega_input"=array(0,dim=c(TmbData$n_knots,TmbData$n_factors)), "delta_i"=rep(0,TmbData$n_obs) )
if(Version=="spatial_dfa_v3") TmbParams = list("logkappa_j"=rep(log(100),TmbData$n_factors), "alpha_j"=rep(0,TmbData$n_factors), "phi_j"=rep(0,TmbData$n_factors), "loglambda_j"=rep(log(1),TmbData$n_factors), "rho_j"=rep(0.2,TmbData$n_factors), "L_val"=rnorm(TmbData$n_factors*TmbData$n_species-TmbData$n_factors*(TmbData$n_factors-1)/2), "gamma_k"=rep(1,TmbData$n_cov), "log_sigma_p"=rep(log(1),TmbData$n_species), "Epsilon_input"=array(0,dim=c(TmbData$n_knots,TmbData$n_factors,TmbData$n_years)), "Omega_input"=array(0,dim=c(TmbData$n_knots,TmbData$n_factors)), "delta_i"=rep(0,TmbData$n_obs) )
if(Version=="spatial_dfa_v4") TmbParams = list("logkappa_j"=rep(log(100),TmbData$n_factors), "alpha_j"=rep(0,TmbData$n_factors), "phi_j"=rep(0,TmbData$n_factors), "loglambda_j"=rep(log(1),TmbData$n_factors), "rho_j"=rep(0.2,TmbData$n_factors), "L_val"=rnorm(TmbData$n_factors*TmbData$n_species-TmbData$n_factors*(TmbData$n_factors-1)/2), "gamma_k"=rep(1,TmbData$n_cov), "log_sigma_p"=rep(log(1),TmbData$n_species), "zinfl_pz"=matrix(0,nrow=TmbData$n_species,ncol=2), "Epsilon_input"=array(0,dim=c(TmbData$n_knots,TmbData$n_factors,TmbData$n_years)), "Omega_input"=array(0,dim=c(TmbData$n_knots,TmbData$n_factors)), "delta_i"=rep(0,TmbData$n_obs) )

# Random                                                                                                                                                                                                                                                                                                                                                                                                                              
Random = c( "Omega_input", "Epsilon_input", "delta_i" )
if(Use_REML==TRUE) Random = c(Random, "gamma_k")  

# Fixed values
Map = list()
Map[["alpha_j"]] = factor( rep(NA,length(TmbParams[["alpha_j"]])) )
# Parameters shared between dynamic factors
Map[["logkappa_j"]] = factor( rep(1,length(TmbParams[["logkappa_j"]])) )
Map[["rho_j"]] = factor( rep(1,length(TmbParams[["rho_j"]])) )
Map[["loglambda_j"]] = factor( rep(1,length(TmbParams[["loglambda_j"]])) )

# Run-specific fixed values
if(ObsModel==0){
  Map[["zinfl_pz"]] = factor( rep(NA,prod(dim(TmbParams[["zinfl_pz"]]))) )
}
if(ObsModel==1){
  Map[["delta_i"]] = factor( rep(NA,length(TmbParams[["delta_i"]])) )
}

# Link TMB 
setwd( TmbFile )
dyn.load( dynlib(Version) )                                                         # log_tau=0.0,

# Initialization
obj <- MakeADFun(data=TmbData, parameters=TmbParams, random=Random, map=Map, hessian=FALSE, inner.control=list(maxit=1000) )
obj$control <- c( obj$control, list(trace=1, parscale=1, REPORT=1, reltol=1e-12, maxit=100) )
obj$env$inner.control <- c(obj$env$inner.control, list("step.tol"=1e-8, "tol10"=1e-6, "grad.tol"=1e-8) )

# First marginal likelihood and gradient
obj$fn( obj$par )
obj$gr( obj$par )

# Bounds
Upper = rep(Inf, length(obj$par) )
  Upper[grep("rho_j",names(obj$par))] = 0.99
Lower = rep(-Inf, length(obj$par) )
  Lower[grep("rho_j",names(obj$par))] = -0.99

# Run model
opt = nlminb(start=obj$env$last.par.best[-c(obj$env$random)], objective=obj$fn, gradient=obj$gr, upper=Upper, lower=Lower, control=list(eval.max=1e4, iter.max=1e4, trace=1, rel.tol=1e-14) )
opt[["final_gradient"]] = obj$gr( opt$par )

# Loadings matrix
L_pj = Report$L_pj
dimnames(L_pj) = list(levels(DF[,'spp']), paste("Factor",1:Nfactors))

# Extract factors
Psi = Report$psi_njt

# Varimax
if(Nfactors>1){
  Hinv = varimax( L_pj, normalize=FALSE )
    L_pj_rot = L_pj %*% Hinv$rotmat
  Psi_rot = array(NA, dim=dim(Psi))
  for( n in 1:Nknots ) Psi_rot[n,,] = solve(Hinv$rotmat) %*% Psi[n,,]
}else{
  L_pj_rot = L_pj
  Psi_rot = Psi
}

# Cross-correlations
CovMat = L_pj_rot %*% t(L_pj_rot) # Identical to: CorrMat = L_pj %*% t(L_pj)
CorrMat = L_pj_rot %*% t(L_pj_rot) / outer( sqrt(rowSums(L_pj_rot^2)), sqrt(rowSums(L_pj_rot^2)))
