
Rotate_Fn = function( L_pj, Psi, RotationMethod="Varimax" ){

  # Varimax
  if( RotationMethod=="Varimax" ){
    Hinv = varimax( L_pj, normalize=FALSE )
      L_pj_rot = L_pj %*% Hinv$rotmat
    Psi_rot = array(NA, dim=dim(Psi))
    for( n in 1:Nknots ) Psi_rot[n,,] = solve(Hinv$rotmat) %*% Psi[n,,]
  }

  # PCA
  if( RotationMethod=="PCA" ){
    Eigen = eigen(L_pj%*%t(L_pj))
    # Check decomposition
    #all(approx_equal( Eigen$vectors%*%diag(Eigen$values)%*%t(Eigen$vectors), L_pj%*%t(L_pj)))
    # My attempt at new loadings matrix
    L_pj_rot = (Eigen$vectors%*%diag(sqrt(trunc_machineprec(Eigen$values))))[,1:Nfactors,drop=FALSE]
    rownames(L_pj_rot) = rownames(L_pj)
    # My new factors
    Hinv = list("rotmat"=corpcor::pseudoinverse(L_pj_rot)%*%L_pj)
    Psi_rot = array(NA, dim=dim(Psi))
    for( n in 1:Nknots ) Psi_rot[n,,] = Hinv$rotmat %*% Psi[n,,]
  }

  # Check for errors
  approx_equal = function(m1,m2,d=0.001) 2*abs(m1-m2)/(m1+m2)<d
  # Check covariance matrix
  if( !all(approx_equal(L_pj%*%t(L_pj),L_pj_rot%*%t(L_pj_rot))) ) stop("Rotation matrix is wrong")
  # Check linear predictor
  for(i in 1:dim(Psi)[[1]]){
  for(j in 1:dim(Psi)[[3]]){
    if( !all(approx_equal(L_pj%*%Psi[i,,j],L_pj_rot%*%Psi_rot[i,,j])) ) stop(paste0("Linear predictor is wrong for site ",i," and time ",j))
  }}
  
  # Return stuff
  Return = list( "L_pj_rot"=L_pj_rot, "Psi_rot"=Psi_rot, "Hinv"=Hinv )
  return( Return )
}
