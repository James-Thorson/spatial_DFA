// Space time 
#include <TMB.hpp>

// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

// dlognorm
template<class Type>
Type dlognorm(Type x, Type meanlog, Type sdlog, int give_log=false){
  Type Return;
  if(give_log==false) Return = dnorm( log(x), meanlog, sdlog, false) / x;
  if(give_log==true) Return = dnorm( log(x), meanlog, sdlog, true) - log(x);
  return Return;
}

// Main function
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Settings
  DATA_FACTOR( Options_vec );
  // Slot 0 -- distribution of data (0=Normal; 1=Lognormal)
  // Slot 1 -- link for linpred prior to summing (0=None; 1=DEPRECATED; 2=log-link)

  // Indices
  DATA_INTEGER(n_obs);       // Total number of observations (i)
  DATA_INTEGER(n_years);	   // Number of years (t)
  DATA_INTEGER(n_species);   // Number of species (p)
  DATA_INTEGER(n_factors);   // Number of dynamic factors (j)

  // Data
  DATA_VECTOR( c_i );         // Count for observation
  DATA_VECTOR( predTF_i );    // Is observation i included in likelihood (0) or predictive score (1)
  DATA_FACTOR( p_i );       	// Species for observation
  DATA_FACTOR( t_i );       	// Year for observation

  // Fixed effects
  PARAMETER_VECTOR(alpha_p);   // Mean of Gompertz-drift field
  PARAMETER_VECTOR(rho_j);             // Autocorrelation (i.e. density dependence)
  PARAMETER_VECTOR(L_val);    // Values in loadings matrix
  PARAMETER_VECTOR(log_sigma_p);

  // Random effects
  PARAMETER_MATRIX(d_jt);

  // Global stuff
  Type pi = 3.141592;
  Type jnll = 0;
  vector<Type> jnll_c(2);
  jnll_c.setZero();

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

  // Computate probability of factors
  for(int j=0; j<n_factors; j++){
    jnll_c(1) -= dnorm( d_jt(j,0), Type(0.0), Type(1.0), true );
    for(int t=1; t<n_years; t++){
      jnll_c(1) -= dnorm( d_jt(j,t), d_jt(j,t-1)*rho_j(j), Type(1.0), true );
    }
  }

  // Computate expected log-densities
  array<Type> theta_pt(n_species, n_years);
  array<Type> linpred_pt(n_species, n_years);
  for(int p=0; p<n_species; p++){
  for(int t=0; t<n_years; t++){
    for(int j=0; j<n_factors; j++){
      if( Options_vec(1)==0 ) theta_pt(p,t) += d_jt(j,t) * L_pj(p,j);
      if( Options_vec(1)==1 ) theta_pt(p,t) += exp( d_jt(j,t) * L_pj(p,j) );          // DEPRECATED
      if( Options_vec(1)==2 ) theta_pt(p,t) += exp( d_jt(j,t) + L_pj(p,j) );
    }
    // Expectation
    if( Options_vec(1)==0 ) linpred_pt(p,t) = theta_pt(p,t) + alpha_p(p);
    if( Options_vec(1)==1 | Options_vec(1)==2 ) linpred_pt(p,t) = log( theta_pt(p,t) + exp(alpha_p(p)) );
  }}

  // Probability of observations
  vector<Type> nll_i(n_obs);
  for(int i=0; i<n_obs; i++){
    // Likelihood
    if( !isNA(c_i(i)) ){                
      if( Options_vec(0)==0 ) nll_i(i) = -1 * dnorm( c_i(i), linpred_pt(p_i(i),t_i(i)), exp(log_sigma_p(p_i(i))), true );
      if( Options_vec(0)==1 & c_i(i)!=0 ) nll_i(i) = -1 * dlognorm( c_i(i), linpred_pt(p_i(i),t_i(i)), exp(log_sigma_p(p_i(i))), true );
    }
  }
  jnll_c(0) = ( (Type(1.0)-predTF_i)*nll_i ).sum();
  Type pred_jnll; // = ( predTF_i*nll_i ).sum();
  jnll = jnll_c.sum();

  // Penalties
  REPORT( nll_i );
  REPORT( jnll_c );
  REPORT( jnll );
  REPORT( pred_jnll );

  // Loadings
  REPORT( L_pj );

  // Predictions
  REPORT( linpred_pt );

  // Derived fields
  REPORT( theta_pt );

  // Parameters
  REPORT( d_jt );
  REPORT( alpha_p );
  REPORT( rho_j );
  REPORT( log_sigma_p );

  // Return objective function
  return jnll;
}
