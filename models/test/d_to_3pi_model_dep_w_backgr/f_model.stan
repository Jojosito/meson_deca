model {
}
generated quantities{
  vector[num_resonances()] A_cv_y[2];
  vector[num_resonances()] theta[2];
  vector[num_background()] A_v_background_abs2_;
  vector[num_background()] theta_background_abs2_;
  real z;

  // Fill the variables
  A_cv_y[1,1] <- 1.0;
  A_cv_y[2,1] <- 0.0;
  A_cv_y[1,2] <- 1.0;
  A_cv_y[2,2] <- 0.0;
  A_cv_y[1,3] <- 1.0;
  A_cv_y[2,3] <- 0.0;
  A_cv_y[1,4] <- 1.0;
  A_cv_y[2,4] <- 0.0;
  A_cv_y[1,5] <- 1.0;
  A_cv_y[2,5] <- 0.0;
  A_cv_y[1,6] <- 1.0;
  A_cv_y[2,6] <- 0.0;
  A_cv_y[1,7] <- 1.0;
  A_cv_y[2,7] <- 0.0;

  theta[1,1] <- 1.0;
  theta[2,1] <- 0.0;
  theta[1,2] <- 1.0;
  theta[2,2] <- 0.0;
  theta[1,3] <- 1.0;
  theta[2,3] <- 0.0;
  theta[1,4] <- 1.0;
  theta[2,4] <- 0.0;
  theta[1,5] <- 1.0;
  theta[2,5] <- 0.0;
  theta[1,6] <- 1.0;
  theta[2,6] <- 0.0;
  theta[1,7] <- 1.0;
  theta[2,7] <- 0.0;

  A_v_background_abs2_[1] <- 1.0;
  A_v_background_abs2_[2] <- 0.0;

  theta_background_abs2_[1] <- 1.0;
  theta_background_abs2_[2] <- 0.0;

  z <- f_model(A_cv_y, theta,
               A_v_background_abs2_, theta_background_abs2_);

  print("f_model(..):");
  print(z);
}