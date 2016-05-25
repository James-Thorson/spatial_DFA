
# Install package
#devtools::install_github("james-thorson/Spatial_DFA")

# File structure
#TmbFile = system.file("executables", package="SpatialDFA")
TmbFile = "C:/Users/James.Thorson/Desktop/Project_git/spatial_DFA/inst/executables"

# Settings
Version = "spatial_dfa_v18"
Sim_Settings = list("n_species"=5, "n_years"=20, "n_stations"=20, "n_factors"=2, "simulation_method"=c("mesh","grid")[2], "SD_O"=0.5, "SD_E"=0.2, "rho"=0.8, "SpatialScale"=NA, "SD_extra"=0.05)
Sim_Settings[["SpatialScale"]] = ifelse( Sim_Settings[["simulation_method"]]=="grid", 2, 0.25 )

# Settings
Nfactors = 2
estimation_method = c("mesh","grid")[2]

# Libraries
library( INLA )
library( TMB )

# Specific libraries
library( SpatialDFA )  

# Compile TMB model
setwd( TmbFile )
compile( paste0(Version,".cpp") )

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
InputList = MakeInput_Fn( Version=Version, DF=DF, Nfactors=Nfactors, loc_xy=loc_xy, method=estimation_method )

# Link TMB 
dyn.load( dynlib(Version) )                                                         # log_tau=0.0,

# Initialization
start_time = Sys.time()
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
for(i in 1:3) opt = nlminb(start=obj$env$last.par.best[-c(obj$env$random)], objective=obj$fn, gradient=obj$gr, upper=Upper, lower=Lower, control=list(eval.max=1e4, iter.max=1e4, trace=1, rel.tol=1e-14) )
opt[["final_gradient"]] = obj$gr( opt$par )
opt[["run_time"]] = Sys.time() - start_time
Report = obj$report()
ParHat = obj$env$parList()

# Loadings matrix
L_pj = Report$L_pj
dimnames(L_pj) = list(levels(DF[,'spp']), paste("Factor",1:Nfactors))

# Extract factors
Psi = Report$psi_njt

# Rotate
if(Nfactors>1){
  RotateList = Rotate_Fn( L_pj=L_pj, Psi=Psi, RotationMethod="PCA", testcutoff=1e-5 )
  L_pj_rot = RotateList[["L_pj_rot"]]
  Psi_rot = RotateList[["Psi_rot"]]
}else{
  L_pj_rot = L_pj
  Psi_rot = Psi
}

##############
# Comparison of true and estimated quantities
##############

# Correlation
cov2cor(L_pj_rot %*% t(L_pj_rot))
cov2cor(Sim_List$Lmat %*% t(Sim_List$Lmat))

# Decorrelation distance
if( Sim_Settings[["simulation_method"]]=="mesh" & estimation_method=="mesh" ){
  print(Sim_Settings[["SpatialScale"]])
  print(Report$Range_jz/2)
}
if( Sim_Settings[["simulation_method"]]=="grid" & estimation_method=="grid" ){
  print(Sim_Settings[["SpatialScale"]])
  print( exp(ParHat$logkappa_jz) )
}

# Ratio of spatial and spatio-temporal variation
(Sim_Settings[["SD_E"]] / Sim_Settings[["SD_O"]])^2
exp(ParHat$loglambda_j)

# Autocorrelation over time
Sim_Settings[["rho"]]
ParHat$rho

