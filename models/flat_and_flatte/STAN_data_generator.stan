functions{
}

data {
  vector[num_resonances()] theta[2];
}

parameters {
  vector<lower=0., upper=3.>[num_variables()] y;
}

model {

  real logH;
  logH <- 0;

  logH <- logH + log( f_model(A_cv(y), theta) );
  increment_log_prob(logH);

}






















