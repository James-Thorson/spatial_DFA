
# Install package
#devtools::install_github("james-thorson/Spatial_DFA")
#devtools::install_github("kaskr/TMB_contrib_R/TMBhelper")

# File structure
TmbFile = system.file("executables", package="SpatialDFA")
#TmbFile = "C:/Users/James.Thorson/Desktop/Project_git/spatial_DFA/inst/executables"

# Libraries
library( INLA )
library( TMB )

# Specific libraries
library( SpatialDFA )
Version = "spatial_dfa_v18"

# Compile TMB model
setwd( TmbFile )
compile( paste0(Version,".cpp") )

# Settings
set.seed(1)
Sim_Settings = list("n_species"=5, "n_years"=20, "n_stations"=49, "n_factors"=3, "simulation_method"=c("mesh","grid")[2], "SD_O"=0.5, "SD_E"=0.2, "rho"=0.8, "SpatialScale"=NA, "SD_extra"=0.05)
Sim_Settings[["SpatialScale"]] = ifelse( Sim_Settings[["simulation_method"]]=="grid", 2, 0.25 )

# Settings
Nfactors = 2
estimation_method = c("mesh","grid")[2]

##############
# Simulation
##############

# Simulate data
Sim_List = Sim_Fn(n_species=Sim_Settings[["n_species"]], simulation_method=Sim_Settings[["simulation_method"]], "SD_O"=Sim_Settings[["SD_O"]], "SD_E"=Sim_Settings[["SD_E"]], "rho"=Sim_Settings[["rho"]], n_years=Sim_Settings[["n_years"]], n_stations=Sim_Settings[["n_stations"]], n_factors=Sim_Settings[["n_factors"]], SpatialScale=Sim_Settings[["SpatialScale"]], SD_extra=Sim_Settings[["SD_extra"]])

# Unload stuff
DF = Sim_List[["DF"]]

# lat_set
loc_xy = Sim_List$Loc

# Bundle inputs
#Nobsfactors=0; Kappa_Type="Constant"; ObsModel=NULL;
#  Aniso=FALSE; Include_Omega=TRUE; Include_Epsilon=TRUE; EncounterFunction=2; Correlated_Overdispersion=FALSE;
#  Include_Phi=TRUE; Include_Rho=TRUE; Use_REML=FALSE; X_ik=NULL; X_nl=NULL; X_ntl=NULL; a_n=NULL; YearSet=NULL;
#  IndependentTF=c(FALSE,FALSE); CheckForBugs=TRUE; CorrGroup_pp=NULL
InputList = MakeInput_Fn( Version=Version, DF=DF, Nfactors=Nfactors, loc_xy=loc_xy, method=estimation_method )

# Link TMB 
dyn.load( dynlib(Version) )

# Initialization
obj <- MakeADFun(data=InputList[["TmbData"]], parameters=InputList[["TmbParams"]], random=InputList[["Random"]], map=InputList[["Map"]], hessian=FALSE, inner.control=list(maxit=1000) )
obj$control <- c( obj$control, list(trace=1, parscale=1, REPORT=1, reltol=1e-12, maxit=100) )
obj$env$inner.control <- c(obj$env$inner.control, list("step.tol"=1e-8, "tol10"=1e-6, "grad.tol"=1e-8) )

# Bounds
Upper = rep(Inf, length(obj$par) )
  Upper[grep("rho_j",names(obj$par))] = 0.99
Lower = rep(-Inf, length(obj$par) )
  Lower[grep("rho_j",names(obj$par))] = -0.99
if( estimation_method=="grid" ){
  Upper[grep("logkappa_jz",names(obj$par))] = -1e-4
}

# Run model
opt = TMBhelper::Optimize(obj=obj, upper=Upper, lower=Lower, getsd=TRUE, newtonsteps=0, loopnum=3 )
Report = obj$report()
ParHat = obj$env$parList()

# Summarize results
Results = Summarize( obj, species_names=levels(DF[,'spp']) )

##############
# Comparison of true and estimated quantities
##############

# Correlation
cov2cor(Results$L_pj_rot %*% t(Results$L_pj_rot))
cov2cor(Sim_List$Lmat %*% t(Sim_List$Lmat))

# Decorrelation distance
print(Sim_Settings[["SpatialScale"]])
print(Report$Range_jz/2)

# Ratio of spatial and spatio-temporal variation
(Sim_Settings[["SD_E"]] / Sim_Settings[["SD_O"]])^2
exp(ParHat$loglambda_j)

# Autocorrelation over time
Sim_Settings[["rho"]]
ParHat$rho

