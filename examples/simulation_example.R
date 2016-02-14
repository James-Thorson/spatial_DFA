
# Install package
devtools::install_github("james-thorson/Spatial_DFA", auth_token=[contact Jim Thorson for a token])

# File structure
TmbFile = system.file("executables", package="SpatialDFA")

# Settings
Version = "spatial_dfa_v17"
Sim_Settings = list("n_species"=5, "n_years"=20, "n_stations"=20, "n_factors"=2, "SpatialScale"=0.25, "SD_extra"=0.05)

# Settings
Nfactors = 2

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
Sim_List = Sim_Fn(n_species=Sim_Settings[["n_species"]], n_years=Sim_Settings[["n_years"]], n_stations=Sim_Settings[["n_stations"]], n_factors=Sim_Settings[["n_factors"]], SpatialScale=Sim_Settings[["SpatialScale"]], SD_extra=Sim_Settings[["SD_extra"]])

# Unload stuff
DF = Sim_List[["DF"]]

# lat_set
latlong_set = paste( Sim_List[["Loc"]][,'y'], Sim_List[["Loc"]][,'x'], sep="_")
lat_set = as.numeric(sapply( latlong_set, FUN=function(Char){strsplit(Char,"_")[[1]][1]}))
long_set = as.numeric(sapply( latlong_set, FUN=function(Char){strsplit(Char,"_")[[1]][2]}))

## Build SPDE object using INLA
mesh = inla.mesh.create( cbind(long_set, lat_set), plot.delay=NULL, extend=list(n=8,offset=-0.15), refine=list(min.angle=26) )  # loc_samp  ;  ,max.edge.data=0.08,max.edge.extra=0.2

# Bundle inputs
InputList = MakeInput_Fn( Version=Version, DF=DF, Nfactors=Nfactors, inla_mesh=mesh )

# Link TMB 
dyn.load( dynlib(Version) )                                                         # log_tau=0.0,

# Initialization
obj <- MakeADFun(data=InputList[["TmbData"]], parameters=InputList[["TmbParams"]], random=InputList[["Random"]], map=InputList[["Map"]], hessian=FALSE, inner.control=list(maxit=1000) )
obj$control <- c( obj$control, list(trace=1, parscale=1, REPORT=1, reltol=1e-12, maxit=100) )
obj$env$inner.control <- c(obj$env$inner.control, list("step.tol"=1e-8, "tol10"=1e-6, "grad.tol"=1e-8) )

# Bounds
Upper = rep(Inf, length(obj$par) )
  Upper[grep("rho_j",names(obj$par))] = 0.99
Lower = rep(-Inf, length(obj$par) )
  Lower[grep("rho_j",names(obj$par))] = -0.99

# Run model
opt = nlminb(start=obj$env$last.par.best[-c(obj$env$random)], objective=obj$fn, gradient=obj$gr, upper=Upper, lower=Lower, control=list(eval.max=1e4, iter.max=1e4, trace=1, rel.tol=1e-14) )
opt[["final_gradient"]] = obj$gr( opt$par )
Report = obj$report()

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

# Cross-correlations
cov2cor(L_pj_rot %*% t(L_pj_rot))

# True cross-correlations
cov2cor(Sim_List$Lmat %*% t(Sim_List$Lmat))
