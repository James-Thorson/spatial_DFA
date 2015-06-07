
InjectParams = function( previous, skeleton ){
  # Start from skeleton
  new = skeleton
  # Identical to previous MLE
  for(i in 1:4 ){
    parname = c("delta_i", "gamma_k", "log_sigma_p", "zinfl_pz")[i]
    new[[parname]] = previous[[parname]]
  }
  # Needs new vector slots
  for(i in 1:5){
    parname = c("alpha_j", "phi_j", "loglambda_j", "rho_j", "L_val")[i]
    new[[parname]][1:length(previous[[parname]])] = previous[[parname]]
  }
  # Needs one new row (identical to the row above because of how "map" input works)
  new[["logkappa_jz"]][1:nrow(previous[["logkappa_jz"]]),][] = as.vector(previous[["logkappa_jz"]])
  new[["logkappa_jz"]][nrow(new[["logkappa_jz"]]),] = previous[["logkappa_jz"]][nrow(previous[["logkappa_jz"]]),]
  # Needs new column
  new[["Omega_input"]][,1:ncol(previous[["Omega_input"]])][1:length(as.vector(previous[["Omega_input"]]))] = as.vector(previous[["Omega_input"]])
  # Needs new column
  new[["Epsilon_input"]][,1:dim(previous[["Epsilon_input"]])[2],][1:length(as.vector(previous[["Epsilon_input"]]))] = as.vector(previous[["Epsilon_input"]])
  # Return new parameter list
  return(new)
}
