
#' Summarize and explain results
#'
#' \code{Summarize} is intended to help document results, by providing explanation for outputs from \code{obj$report()}
#'
#' @param obj the TMB object after it has been previously optimized
#' @param species_names an optional character-vector giving species names
#' @param rotation_method what rotation method is used when visualizing results (default \code{rotation_method="PCA"})

#' @return Tagged list of useful output
#' \describe{
#'   \item{parhat}{a list of fixed-effect estimates, with same format as starting values}
#'   \item{L_pj}{the loadings matrix, where e.g. \code{L_pj[1,2]} gives the association of species 1 with factor 2}
#'   \item{Psi}{the array of factors, where e.g. \code{Psi[1,2,3]} gives the value at site 1 of factor 2 in time 3}
#'   \item{L_pj_rot}{the loadings matrix after rotation, useful for ease of interpretation}
#'   \item{Psi_rot}{the value of factors, useful for ease of interpretation}
#'   \item{range_jz}{Decorrelation distance for factor j (1st column: spatial; 2nd column: spatio-temporal}
#'   \item{theta_npt}{Empirical Bayes prediction of log-density, where e.g., \code{theta_npt[1,2,3]} gives the log-density for site 1 and species 2 in time 3}
#'   \item{eta_npt}{Empirical Bayes prediction of log-detectability, where e.g., \code{eta_npt[1,2,3]} gives the log-detectability for site 1 and species 2 in time 3}
#'   \item{logchat_i}{Empirical Bayes prediction of log-expectation for observation i}
#'   \item{nll_i}{Empirical Bayes prediction of one-half of the deviance residual for observation i}
#' }

#' @export
Summarize = function( obj, species_names=NULL, rotation_method="PCA" ){
  # fill in defaults
  if( is.null(species_names) ) species_names = paste0("Species_",1:obj$env$data$n_species)

  # extract elements
  report = obj$report()
  parhat = obj$env$parList()

  # Loadings matrix
  L_pj = report$L_pj
  dimnames(L_pj) = list(species_names, paste0("Factor_",1:obj$env$data$n_factors))

  # Extract factors
  Psi = report$psi_njt

  # Rotate
  if(Nfactors>1){
    RotateList = SpatialDFA:::Rotate_Fn( L_pj=L_pj, Psi=Psi, RotationMethod=rotation_method, testcutoff=1e-5 )
    L_pj_rot = RotateList[["L_pj_rot"]]
    Psi_rot = RotateList[["Psi_rot"]]
  }else{
    L_pj_rot = L_pj
    Psi_rot = Psi
  }

  # Bundle and return
  Return = list( "parhat"=parhat, "L_pj"=L_pj, "Psi"=Psi, "L_pj_rot"=L_pj_rot, "Psi_rot"=Psi_rot, "range_jz"=report$Range_jz, "theta_npt"=report$theta_npt, "eta_npt"=report$eta_npt, "logchat_i"=report$logchat_i, "nll_i"=report$nll_i)
  return( Return )
}

