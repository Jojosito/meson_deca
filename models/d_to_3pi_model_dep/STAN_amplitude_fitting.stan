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
  // Total: 12
  real<lower=-5., upper=5.> theta_re_flat;
  real<lower=-5., upper=5.> theta_im_flat;

  real<lower=-5., upper=5.> theta_re_f0_980;
  real<lower=-5., upper=5.> theta_im_f0_980;

  real<lower=-5., upper=5.> theta_re_f0_600;
  real<lower=-5., upper=5.> theta_im_f0_600;

  real<lower=-5., upper=5.> theta_re_f0_1370;
  real<lower=-5., upper=5.> theta_im_f0_1370;

  real<lower=-5., upper=5.> theta_re_f0_1500;
  real<lower=-5., upper=5.> theta_im_f0_1500;

  real<lower=-5., upper=5.> theta_re_f2_1270;
  real<lower=-5., upper=5.> theta_im_f2_1270;


}

transformed parameters {

  // Parameters: some fixed (reference parameters), some free
  // (these will be fitted).
  vector<lower=-5., upper=5.>[num_resonances()] theta[2];

  theta[1,1] <- theta_re_flat;
  theta[2,1] <- theta_im_flat;

  theta[1,2] <- theta_re_f0_980;
  theta[2,2] <- theta_im_f0_980;

  theta[1,3] <- theta_re_f0_600;
  theta[2,3] <- theta_im_f0_600;

  theta[1,4] <- theta_re_f0_1370;
  theta[2,4] <- theta_im_f0_1370;

  theta[1,5] <- theta_re_f0_1500;
  theta[2,5] <- theta_im_f0_1500;

  theta[1,6] <- 1.0;
  theta[2,6] <- 0.0;

  theta[1,7] <- theta_re_f2_1270;
  theta[2,7] <- theta_im_f2_1270;

}


model {

  real logH;
  logH <- 0;

  // Sum over all events
  for (d in 1:D)
      logH <- logH + log( f_model(A_cv_data[d], theta) / Norm(theta, I) );

  increment_log_prob(logH);

}
