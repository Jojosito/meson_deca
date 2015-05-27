functions{
}


data {

  // Number of measured events
  int D;

  // Complex PWA amplitudes corresponding to each event
  vector[num_resonances()] A_cv_data[D,2];

  // Complex normalization matrix corresponding to the model
  matrix[num_resonances(), num_resonances()] I[2];

}


parameters {
  // Parameters that will be fitted
  real<lower=-2., upper=2.> theta_re;
  real<lower=-2., upper=2.> theta_im;
}

transformed parameters {

  // Parameters: some fixed (reference parameters), some free
  // (these will be fitted).
  vector<lower=-2., upper=2.>[num_resonances()] theta[2];

  theta[1,1] <- 1.0;
  theta[2,1] <- 0.0;
  theta[1,2] <- theta_re;
  theta[2,2] <- theta_im;
  theta[1,3] <- 0.0;
  theta[2,3] <- 0.0;

}


model {

  real logH;
  logH <- 0;

  // Sum over all events
  for (d in 1:D)
      logH <- logH + log( f_model(A_cv_data[d], theta) / Norm(theta, I) );

  increment_log_prob(logH);

}
