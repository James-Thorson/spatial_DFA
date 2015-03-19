// Space time 
#include <TMB.hpp>
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}
template<class Type>
Type objective_function<Type>::operator() ()
{
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

  // SPDE objects
  DATA_SPARSE_MATRIX(G0);
  DATA_SPARSE_MATRIX(G1);
  DATA_SPARSE_MATRIX(G2);

  // Fixed effects
  PARAMETER_VECTOR(logkappa_j);         // Controls range of spatial variation
  PARAMETER_VECTOR(alpha_j);   // Mean of Gompertz-drift field
  PARAMETER_VECTOR(phi_j);              // Offset of beginning from equilibrium
  PARAMETER_VECTOR(loglambda_j);        // ratio of temporal and spatial variance for dynamic-factor j
  PARAMETER_VECTOR(rho_j);             // Autocorrelation (i.e. density dependence)
  PARAMETER_VECTOR(L_val);    // Values in loadings matrix
  PARAMETER_VECTOR(gamma_k);
  PARAMETER_VECTOR(log_sigma_p);

  // Random effects
  PARAMETER_ARRAY(Epsilon_input);  // Spatial process variation
  PARAMETER_MATRIX(Omega_input);   // Spatial variation in carrying capacity
  PARAMETER_VECTOR(delta_i);

  // Global stuff
  using namespace density;
  Type pi = 3.141592;
  Type jnll = 0;

  // Calculate tau (space: O, time: E) for dynamics factors
  vector<Type> lambda_j(n_factors);
  vector<Type> kappa_j(n_factors);
  vector<Type> logtau_E_j(n_factors);        // log-inverse SD of Epsilon
  vector<Type> logtau_O_j(n_factors);        // log-inverse SD of Omega
  vector<Type> VarSpace_j(n_factors);
  vector<Type> VarTime_j(n_factors);
  vector<Type> kappa2_j(n_factors);
  vector<Type> kappa4_j(n_factors);
  vector<Type> Range_j(n_factors);
  for(int j=0; j<n_factors; j++){
    kappa2_j(j) = exp(2.0*logkappa_j(j));
    kappa4_j(j) = kappa2_j(j)*kappa2_j(j);
    Range_j(j) = sqrt(8.0) / exp( logkappa_j(j) );
    lambda_j(j) = exp(loglambda_j(j));
    kappa_j(j) = exp(logkappa_j(j));
    logtau_O_j(j) = 0.5*( log(1.0+lambda_j(j)) - log(4.0*pi*pow(kappa_j(j),2.0)) );
    logtau_E_j(j) = logtau_O_j(j) - 0.5*log( lambda_j(j) * (1-pow(rho_j(j),2.0)) );
    VarSpace_j(j) = 1.0 / ( 4.0*pi*exp(2.0*logtau_O_j(j))*exp(2.0*logkappa_j(j)) );
    VarTime_j(j) = 1.0 / ( 4.0*pi*exp(2.0*logtau_E_j(j))*exp(2.0*logkappa_j(j)) * (1-pow(rho_j(j),2.0)) );
  }

  // Assemble the loadings matrix
  matrix<Type> L_pj(n_species,n_factors);
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
  
  // Derived quantities related to GMRF
  Eigen::SparseMatrix<Type> Q;

  // Probability of random fields
  array<Type> Epsilon_tmp(n_knots,n_years);
  vector<Type> pen_epsilon_j(n_factors);
  vector<Type> pen_omega_j(n_factors);
  for(int j=0; j<n_factors; j++){
    Q = kappa4_j(j)*G0 + Type(2.0)*kappa2_j(j)*G1 + G2;
    for(int n=0; n<n_knots; n++){
    for(int t=0; t<n_years; t++){
      Epsilon_tmp(n,t) = Epsilon_input(n,j,t);
    }}
    pen_epsilon_j(j) += SEPARABLE(AR1(rho_j(j)),GMRF(Q))(Epsilon_tmp);
    pen_omega_j(j) += GMRF(Q)(Omega_input.col(j));
  }
  jnll += sum( pen_epsilon_j );
  jnll += sum( pen_omega_j );

  // Probability of overdispersion
  vector<Type> prob_delta_i(n_obs);
  for(int i=0; i<n_obs; i++){
    prob_delta_i(i) = dnorm( delta_i(i), Type(0.0), Type(1.0), true );
  }
  jnll -= sum( prob_delta_i );

  // Transform random fields
  array<Type> Epsilon_njt(n_knots,n_factors,n_years);
  array<Type> Omega_nj(n_knots,n_factors);
  for(int n=0; n<n_knots; n++){
  for(int j=0; j<n_factors; j++){
    Omega_nj(n,j) = Omega_input(n,j) / exp(logtau_O_j(j));
    for(int t=0; t<n_years; t++){
      Epsilon_njt(n,j,t) = Epsilon_input(n,j,t) / exp(logtau_E_j(j));
    }
  }}
  
  // Declare derived quantities
  vector<Type> logchat_i(n_obs);
  vector<Type> mu_i(n_obs);
  array<Type> psi_njt(n_knots,n_factors,n_years); 
  vector<Type> eta_i(n_obs);
  array<Type> theta_npt(n_knots,n_species,n_years);
  
  // Computate factors
  for(int n=0; n<n_knots; n++){
  for(int j=0; j<n_factors; j++){
  for(int t=0; t<n_years; t++){
    psi_njt(n,j,t) = phi_j(j)*pow(rho_j(j),t) + Epsilon_njt(n,j,t) + ( alpha_j(j) + Omega_nj(n,j) )/(1-rho_j(j));
  }}}
  eta_i = X_ik * gamma_k.matrix();

  // Computate expected log-densities
  for(int n=0; n<n_knots; n++){
  for(int p=0; p<n_species; p++){
  for(int t=0; t<n_years; t++){
    theta_npt(n,p,t) = 0.0;
    for(int j=0; j<n_factors; j++){
      theta_npt(n,p,t) += psi_njt(n,j,t) * L_pj(p,j);
    }
  }}}

  // Probability of observations
  for(int i=0; i<n_obs; i++){
    // Expectation
    mu_i(i) = eta_i(i) + theta_npt(s_i(i),p_i(i),t_i(i));
    // Overdispersion
    logchat_i(i) = mu_i(i) + delta_i(i) * exp(log_sigma_p(p_i(i)));
    // Likelihood
    if( !isNA(c_i(i)) ){                
      jnll -= dpois( c_i(i), exp(logchat_i(i)), true );
    }
  }

  // Spatial field summaries
  REPORT( Range_j );
  REPORT( logtau_O_j );
  REPORT( logtau_E_j );
  REPORT( VarSpace_j );
  REPORT( VarTime_j );
  REPORT( rho_j );
  REPORT( logkappa_j );
  REPORT( alpha_j );
  REPORT( phi_j );
  REPORT( lambda_j );
  // Penalties
  REPORT( pen_epsilon_j );
  REPORT( pen_omega_j );
  REPORT( prob_delta_i );
  // Loadings
  REPORT( L_pj );
  // Fields
  REPORT( Epsilon_njt );
  REPORT( Omega_nj );
  // Predictions 
  REPORT( logchat_i );
  REPORT( eta_i );
  REPORT( mu_i );
  // Derived fields
  REPORT( theta_npt );
  REPORT( psi_njt );
  // Parameters
  REPORT( logkappa_j );
  REPORT( alpha_j );
  REPORT( phi_j );
  REPORT( loglambda_j );
  REPORT( rho_j );
  REPORT( gamma_k );
  REPORT( log_sigma_p );
  
  // Return objective function
  return jnll;
}
