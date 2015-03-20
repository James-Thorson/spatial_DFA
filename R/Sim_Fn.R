Sim_Fn <-
function( n_species, n_years, n_stations=20, phi=NULL, n_factors=2, SpatialScale=0.1, SD_O=0.5, SD_E=0.2, SD_extra=0.1, rho=0.8, logMeanDens=1, Lmat=NULL, Loc=NULL ){
  # Parameters
  if( is.null(Lmat) ){
    Lmat = matrix( rnorm(n_factors*n_species), nrow=n_species, ncol=n_factors)
    for(i in 1:ncol(Lmat)){
      Lmat[seq(from=1,to=i-1,length=i-1),i] = 0
      if( Lmat[,i][which.max(abs(Lmat[,i]))]<0 ){
        Lmat[,i] = -1*Lmat[,i]
      }
    }
  }
  if( is.null(phi) ) phi = rnorm(n_factors, mean=0, sd=1)
  Beta = rep(logMeanDens, n_species)

  # Spatial model
  if( is.null(Loc) ) Loc = cbind( "x"=runif(n_stations, min=0,max=1), "y"=runif(n_stations, min=0,max=1) )
  model_O <- RMgauss(var=SD_O^2, scale=SpatialScale)
  model_E <- RMgauss(var=SD_E^2, scale=SpatialScale)

  # Simulate Omega
  Omega = matrix(NA, nrow=n_stations, ncol=n_factors)
  for(i in 1:n_factors){
    Omega[,i] = RFsimulate(model = model_O, x=Loc[,'x'], y=Loc[,'y'])@data[,1]
  }

  # Simulate Epsilon
  Epsilon = array(NA, dim=c(n_stations,n_factors,n_years))
  for(i in 1:n_factors){
    Epsilon[,i,1] = RFsimulate(model=model_E, x=Loc[,'x'], y=Loc[,'y'])@data[,1]
    for(t in 2:n_years){
      Epsilon[,i,t] = rho * Epsilon[,i,t-1] + RFsimulate(model=model_E, x=Loc[,'x'], y=Loc[,'y'])@data[,1]
    }
  }
  
  # Calculate Psi
  Psi = array(NA, dim=c(n_stations,n_factors,n_years))
  for(i in 1:n_factors){
  for(t in 1:n_years){
    Psi[,i,t] = phi[i] * rho^t + Epsilon[,i,t] + Omega[,i]/(1-rho)
  }}
  
  # Expected log-densities
  Theta = array(NA, dim=c(n_stations,n_species,n_years))
  for(s in 1:n_stations){
  for(t in 1:n_years){
    Theta[s,,t] = Lmat %*% Psi[s,,t]
  }}
  
  # Simulate data
  DF = NULL
  for(s in 1:n_stations){
  for(p in 1:n_species){
  for(t in 1:n_years){
    Tmp = c("sitenum"=s, "spp"=p, "year"=t, "catch"=rpois(1,lambda=exp(Theta[s,p,t]+logMeanDens+SD_extra*rnorm(1))), 'waterTmpC'=0 )
    DF = rbind(DF, Tmp)
  }}}
  DF = data.frame(DF, row.names=NULL)
  DF[,'spp'] = factor( letters[DF[,'spp']] )
  if( n_species>26 ) stop( "problem with using letters")

  # Return stuff
  Sim_List = list("DF"=DF, "Psi"=Psi, "Lmat"=Lmat, "phi"=phi, "Loc"=Loc, "Omega"=Omega, "Epsilon"=Epsilon, "Theta"=Theta, "Psi"=Psi)
  return(Sim_List)
}
