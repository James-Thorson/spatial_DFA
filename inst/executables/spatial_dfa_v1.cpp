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
  PARAMETER(log_kappa);         // Controls range of spatial variation
  PARAMETER_VECTOR(alpha_j);   // Mean of Gompertz-drift field
  PARAMETER_VECTOR(phi_j);              // Offset of beginning from equilibrium
  PARAMETER_VECTOR(logtau_E_j);        // log-inverse SD of Epsilon
  PARAMETER_VECTOR(logtau_O_j);        // log-inverse SD of Omega
  PARAMETER_VECTOR(rho_j);             // Autocorrelation (i.e. density dependence)
  PARAMETER_VECTOR(Psi_val);    // Values in loadings matrix
  PARAMETER_VECTOR(gamma_k);
  PARAMETER_VECTOR(log_sigma_p);

  // Random effects
  PARAMETER_ARRAY(Epsilon_input);  // Spatial process variation
  PARAMETER_MATRIX(Omega_input);   // Spatial variation in carrying capacity
  PARAMETER_VECTOR(tilda_i);

  // Assemble the loadings matrix
  matrix<Type> Psi_jp(n_factors,n_species);
  int Count = 0;
  for(int j=0; j<n_factors; j++){
    for(int p=0; p<n_species; p++){
      if(p>=j){
        Psi_jp(j,p) = Psi_val(Count);
        Count++;
      }else{
        Psi_jp(j,p) = 0.0;
      }
    }
  }
  
  // global stuff
  using namespace density;
  Type jnll = 0;
  
  // Derived quantities related to GMRF
  Type kappa2 = exp(2.0*log_kappa);
  Type kappa4 = kappa2*kappa2;
  Type pi = 3.141592;
  Type Range = sqrt(8) / exp( log_kappa );
  Eigen::SparseMatrix<Type> Q = kappa4*G0 + Type(2.0)*kappa2*G1 + G2;

  // Probability of random fields
  array<Type> Epsilon_tmp(n_knots,n_years);
  vector<Type> pen_epsilon_j(n_factors);
  vector<Type> pen_omega_j(n_factors);
  for(int j=0; j<n_factors; j++){
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
  vector<Type> prob_tilda_i(n_obs);
  for(int i=0; i<n_obs; i++){
    prob_tilda_i(i) = dnorm( tilda_i(i), Type(0.0), Type(1.0), true );
  }
  jnll -= sum( prob_tilda_i );

  // Transform random fields
  array<Type> Epsilon(n_knots,n_factors,n_years);
  array<Type> Omega(n_knots,n_factors);
  for(int n=0; n<n_knots; n++){
  for(int j=0; j<n_factors; j++){
    Omega(n,j) = Omega_input(n,j) / exp(logtau_O_j(j));
    for(int t=0; t<n_years; t++){
      Epsilon(n,j,t) = Epsilon_input(n,j,t) / exp(logtau_E_j(j));
    }
  }}
  
  // Declare derived quantities
  vector<Type> logchat_i(n_obs);
  array<Type> lambda_sjt(n_knots,n_factors,n_years); 
  vector<Type> eta_i(n_obs);
  
  // Computate derived quantities
  for(int s=0; s<n_sites; s++){
  for(int j=0; j<n_factors; j++){
  for(int t=0; t<n_years; t++){
    lambda_sjt(s,j,t) = phi_j(j)*pow(rho_j(j),t) + Epsilon(s,j,t) + ( alpha_j(j) + Omega(s,j) )/(1-rho_j(j));
  }}}
  eta_i = X_ik * gamma_k.matrix();

  // Probability of observations
  for(int i=0; i<n_obs; i++){
    logchat_i(i) = eta_i(i);
    for(int j=0; j<n_factors; j++){
      logchat_i(i) += lambda_sjt(s_i(i),j,t_i(i)) * Psi_jp(j,p_i(i));
    }
    // Overdispersion
    logchat_i(i) += tilda_i(i) * exp(log_sigma_p(p_i(i)));
    if( !isNA(c_i(i)) ){                
      jnll -= dpois( c_i(i), exp(logchat_i(i)), true );
    }
  }

  // Spatial field summaries
  REPORT( Range );
  // Fields
  REPORT( Epsilon );
  REPORT( Omega );
  REPORT( lambda_sjt );
  REPORT( logchat_i );
  REPORT( eta_i );
  REPORT( pen_epsilon_j );
  REPORT( pen_omega_j );
  REPORT( Psi_jp );
  REPORT( prob_tilda_i );
  
  return jnll;
}
