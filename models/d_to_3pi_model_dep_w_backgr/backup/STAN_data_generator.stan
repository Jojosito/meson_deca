functions{
}

data {
  // Vector of complex (coherently summed) amplitudes
  vector[num_resonances()] theta[2];
  // Vector of real (incoherently summed) amplitudes
  vector[num_background()] theta_background_abs2;
}

parameters {
  vector<lower=0., upper=3.>[num_variables()] y;
}

model {

  real logH;
  logH <- 0;

  logH <- logH + log( f_model(A_cv(y), theta, 
       	       	      	      A_v_background_abs2(y), theta_background_abs2) );
  increment_log_prob(logH);

}






















