// Space time 
#include <TMB.hpp>

// 2nd power of a number
template<class Type>
Type square(Type x){ return pow(x,2.0); }

// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

// function for logistic transform
template<class Type>
Type plogis(Type x){
  return 1.0 / (1.0 + exp(-x));
}

// dlognorm
template<class Type>
Type dlognorm(Type x, Type meanlog, Type sdlog, int give_log=false){
  Type Return;
  if(give_log==false) Return = dnorm( log(x), meanlog, sdlog, false) / x;
  if(give_log==true) Return = dnorm( log(x), meanlog, sdlog, true) - log(x);
  return Return;
}

// dzinflognorm
template<class Type>
Type dzinflognorm(Type x, Type meanlog, Type sdlog, Type zinfl1, Type zinfl2, int give_log=0){
  Type Return;
  Type encounter_prob = plogis( zinfl1 + zinfl2*meanlog );
  if(x==0){
    if(give_log==false) Return = 1.0 - encounter_prob;
    if(give_log==true) Return = log(1.0 - encounter_prob);
  }else{
    if(give_log==false) Return = encounter_prob * dlognorm( x, meanlog-log(encounter_prob), sdlog, false );
    if(give_log==true) Return = log(encounter_prob) + dlognorm( x, meanlog-log(encounter_prob), sdlog, true );
  } 
  return Return;
}

// Main function
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Settings
  DATA_FACTOR( Options_vec );
  // Slot 0 -- distribution of data
  // Slot 1 -- Include Omega? (0=No, 1=Yes)
  // Slot 2 -- Include Epsilon? (0=No, 1=Yes)
  
  // Indices
  DATA_INTEGER(n_obs);       // Total number of observations (i)
  DATA_INTEGER(n_sites);	   // Number of stations (s)
  DATA_INTEGER(n_years);	   // Number of years (t)
  DATA_INTEGER(n_knots);	   // Number of stations (n)
  DATA_INTEGER(n_species);   // Number of species (p)
  DATA_INTEGER(n_factors);   // Number of dynamic factors (j)
  DATA_INTEGER(n_cov);       // Number of covariates (k)

  // Data
  DATA_VECTOR( c_i );         // Count for observation
  DATA_FACTOR( p_i );       	// Species for observation
  DATA_FACTOR( s_i );       	// Site for observation
  DATA_FACTOR( t_i );       	// Year for observation
  DATA_MATRIX( X_ik );		    // Covariate design matrix

  // Rotation matrix
  DATA_MATRIX( Rotation_jj );

  // SPDE objects
  DATA_SPARSE_MATRIX(G0);
  DATA_SPARSE_MATRIX(G1);
  DATA_SPARSE_MATRIX(G2);

  // Fixed effects
  PARAMETER_MATRIX(logkappa_jz);         // Controls range of spatial variation (0=Omega, 1=Epsilon)
  PARAMETER_VECTOR(alpha_j);   // Mean of Gompertz-drift field
  PARAMETER_VECTOR(phi_j);              // Offset of beginning from equilibrium
  PARAMETER_VECTOR(loglambda_j);        // ratio of temporal and spatial variance for dynamic-factor j
  PARAMETER_VECTOR(rho_j);             // Autocorrelation (i.e. density dependence)
  PARAMETER_VECTOR(L_val);    // Values in loadings matrix
  PARAMETER_VECTOR(gamma_k);
  PARAMETER_VECTOR(log_sigma_p);
  PARAMETER_MATRIX(zinfl_pz);

  // Random effects
  PARAMETER_ARRAY(Epsilon_input);  // Spatial process variation
  PARAMETER_MATRIX(Omega_input);   // Spatial variation in carrying capacity
  PARAMETER_VECTOR(delta_i);

  // Global stuff
  using namespace density;
  Type pi = 3.141592;
  Type jnll = 0;
  vector<Type> jnll_c(4);
  jnll_c.setZero();

  // Calculate tau (space: O, time: E) for dynamics factors
  vector<Type> lambda_j(n_factors);
  vector<Type> VarSpace_j(n_factors);
  vector<Type> VarTime_j(n_factors);
  array<Type> logtau_jz(n_factors,2);        // log-inverse SD of Epsilon
  array<Type> Range_jz(n_factors,2);
  for(int j=0; j<n_factors; j++){
    lambda_j(j) = exp( loglambda_j(j) );
    for(int z=0; z<2; z++) Range_jz(j,z) = sqrt(8.0) / exp( logkappa_jz(j,z) );
    logtau_jz(j,0) = 0.5*( log(1.0+lambda_j(j)) - log(4.0*pi*exp(2.0*logkappa_jz(j,0))*square(1-rho_j(j))) );
    logtau_jz(j,1) = logtau_jz(j,0) + log(1-rho_j(j)) + logkappa_jz(j,0) - logkappa_jz(j,1) - 0.5*log( lambda_j(j) * (1-square(rho_j(j))) );
    VarSpace_j(j) = 1.0 / ( 4.0*pi*exp(2.0*logtau_jz(j,0))*exp(2.0*logkappa_jz(j,0)) * square(1-rho_j(j)) );
    VarTime_j(j) = 1.0 / ( 4.0*pi*exp(2.0*logtau_jz(j,1))*exp(2.0*logkappa_jz(j,1)) * (1-square(rho_j(j))) );
  }

  // Assemble the loadings matrix
  matrix<Type> L_pj(n_species,n_factors);
  matrix<Type> L_pj_rotated(n_species,n_factors);
  int Count = 0;
  for(int j=0; j<n_factors; j++){
  for(int p=0; p<n_species; p++){
    if(p>=j){
      L_pj(p,j) = L_val(Count);
      Count++;
    }else{
      L_pj(p,j) = 0.0;
    }
  }}
  L_pj_rotated = L_pj * Rotation_jj;
  
  // Derived quantities related to GMRF
  Eigen::SparseMatrix<Type> Q;
  // Spatial field reporting
  //GMRF_t<Type> nll_gmrf_spatial(Q)  

  // Probability of random fields
  array<Type> Epsilon_tmp(n_knots,n_years);
  vector<Type> pen_epsilon_j(n_factors);
  vector<Type> pen_omega_j(n_factors);
  for(int j=0; j<n_factors; j++){
    // Omega
    Q = exp(4.0*logkappa_jz(j,0))*G0 + 2.0*exp(2.0*logkappa_jz(j,0))*G1 + G2;
    pen_omega_j(j) += GMRF(Q)(Omega_input.col(j));
    // Epsilon 
    for(int n=0; n<n_knots; n++){
    for(int t=0; t<n_years; t++){
      Epsilon_tmp(n,t) = Epsilon_input(n,j,t);
    }}
    Q = exp(4.0*logkappa_jz(j,1))*G0 + 2.0*exp(2.0*logkappa_jz(j,1))*G1 + G2;
    pen_epsilon_j(j) += SEPARABLE(AR1(rho_j(j)),GMRF(Q))(Epsilon_tmp);
  }
  if( Options_vec(1)==1 ) jnll_c(0) = pen_omega_j.sum();
  if( Options_vec(2)==1 ) jnll_c(1) = pen_epsilon_j.sum();

  // Probability of overdispersion
  vector<Type> prob_delta_i(n_obs);
  for(int i=0; i<n_obs; i++){
    prob_delta_i(i) = dnorm( delta_i(i), Type(0.0), Type(1.0), true );
  }
  jnll_c(2) = -1 * prob_delta_i.sum();

  // Transform random fields
  array<Type> Epsilon_njt(n_knots,n_factors,n_years);
  array<Type> Omega_nj(n_knots,n_factors);
  for(int n=0; n<n_knots; n++){
  for(int j=0; j<n_factors; j++){
    Omega_nj(n,j) = Omega_input(n,j) / exp(logtau_jz(j,0));
    for(int t=0; t<n_years; t++){
      Epsilon_njt(n,j,t) = Epsilon_input(n,j,t) / exp(logtau_jz(j,1));
    }
  }}
  
  // Declare derived quantities
  vector<Type> logchat_i(n_obs);
  array<Type> psi_njt(n_knots,n_factors,n_years); 
  vector<Type> eta_i(n_obs);
  array<Type> theta_npt(n_knots,n_species,n_years);
  
  // Computate factors
  for(int n=0; n<n_knots; n++){
  for(int j=0; j<n_factors; j++){
  for(int t=0; t<n_years; t++){
    psi_njt(n,j,t) = phi_j(j)*pow(rho_j(j),t) + Epsilon_njt(n,j,t) + alpha_j(j)/(1-rho_j(j)) + Omega_nj(n,j)/(1-rho_j(j));
  }}}
  eta_i = X_ik * gamma_k.matrix();

  // Computate expected log-densities
  theta_npt.setZero();
  for(int n=0; n<n_knots; n++){
  for(int p=0; p<n_species; p++){
  for(int t=0; t<n_years; t++){
    for(int j=0; j<n_factors; j++){
      theta_npt(n,p,t) += psi_njt(n,j,t) * L_pj(p,j);
    }
  }}}

  // Probability of observations
  vector<Type> nll_i(n_obs);
  for(int i=0; i<n_obs; i++){
    // Expectation
    logchat_i(i) = eta_i(i) + theta_npt(s_i(i),p_i(i),t_i(i));
    // Overdispersion
    if( Options_vec(0)==0 ) logchat_i(i) += delta_i(i) * exp(log_sigma_p(p_i(i)));
    // Likelihood
    if( !isNA(c_i(i)) ){                
      if( Options_vec(0)==0 ) nll_i(i) = dpois( c_i(i), exp(logchat_i(i)), true );
      if( Options_vec(0)==1 ) nll_i(i) = dzinflognorm( c_i(i), logchat_i(i), exp(log_sigma_p(p_i(i))), zinfl_pz(p_i(i),0), zinfl_pz(p_i(i),1), true );
    }
  }
  jnll_c(3) = -1 * nll_i.sum();
  jnll = jnll_c.sum();

  // Spatial field summaries
  REPORT( Range_jz );
  REPORT( logtau_jz );
  REPORT( VarSpace_j );
  REPORT( VarTime_j );
  REPORT( lambda_j );
  REPORT( Q );
  // Penalties
  REPORT( pen_epsilon_j );
  REPORT( pen_omega_j );
  REPORT( prob_delta_i );
  REPORT( nll_i );
  REPORT( jnll_c );
  REPORT( jnll );
  // Loadings
  REPORT( L_pj );
  REPORT( L_pj_rotated );
  // Fields
  REPORT( Epsilon_njt );
  REPORT( Omega_nj );
  REPORT( Epsilon_input );
  REPORT( Omega_input );
  // Predictions 
  REPORT( logchat_i );
  REPORT( eta_i );
  REPORT( delta_i );
  // Derived fields
  REPORT( theta_npt );
  REPORT( psi_njt );
  // Parameters
  REPORT( logkappa_jz );
  REPORT( alpha_j );
  REPORT( phi_j );
  REPORT( loglambda_j );
  REPORT( rho_j );
  REPORT( gamma_k );
  REPORT( log_sigma_p );
  REPORT( zinfl_pz );
  
  // SEs
  ADREPORT( L_pj_rotated ); 
  
  // Return objective function
  return jnll;
}
