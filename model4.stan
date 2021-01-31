data {
  int<lower=0> N;//Amount of observations
  int observed[N];//Whether the bird was observed or not
  int<lower=0> N_habitats;
  int<lower=0> N_localities;
  vector[N_habitats] habitat_prop[N_localities];
  int localities[N];//The ids have to be preprocessed go 1..N_localities
  vector[N] duration_minutes;
  vector[N] distance_traveled;//Perhaps should also have the protocol type?
  vector[N] forest_prop;//Used for detection model
  vector[N] start_time;
  int<lower=0> check_N;
  int check_indices_habprop[check_N];
  int check_indices_others[check_N];
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  //Obs. model parameters
  vector[N_habitats] habitat_coef;
  real habitat_bias;
  
  //Detection model parameters
  real duration_coef;
  real distance_coef;
  real forest_prop_coef;
  real detectability_bias;
  //Start time model parameters
  real a1;
  real b1;
  real a2;
  real b2;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  vector[N_localities] occupied_prob;
  vector[N] detected_prob;
  
  //Occupancy model
  for(j in 1:N_habitats){//priors for coefs
    habitat_coef[j] ~ normal(0,15);
  }
  habitat_bias ~ normal(0,6);
  for(i in 1:N_localities){
    occupied_prob[i] = inv_logit(habitat_bias + sum(habitat_coef .* habitat_prop[i,]));
  }
  
  //Detection model
  duration_coef ~ normal(0,1.5);
  distance_coef ~ normal(0,6);
  forest_prop_coef ~ normal(0,15);
  detectability_bias ~ normal(0,6);
  a1 ~ normal(0,6);
  b1 ~ normal(0,6);
  a2 ~ normal(0,6);
  b2 ~ normal(0,6);
  for(k in 1:N){
    detected_prob[k] = inv_logit(detectability_bias + duration_coef*duration_minutes[k] + 
                                distance_coef*distance_traveled[k] + 
                                forest_prop_coef*forest_prop[k] + 
                                a1*cos(2*pi()*start_time[k]/24) + b1*sin(2*pi()*start_time[k]/24) + 
                                a2*cos(2*pi()*2*start_time[k]/24) + b2*sin(2*pi()*2*start_time[k]/24));
  }
  
  //Combining the two
  for(k in 1:N){
    observed[k] ~ bernoulli(occupied_prob[localities[k]] * detected_prob[k]);
  }
}

generated quantities{
  int observed_check[check_N];
  vector[check_N] occupied_prob;
  vector[check_N] detected_prob;
  
  for(i in 1:check_N){
    occupied_prob[i] = inv_logit(habitat_bias + sum(habitat_coef .* habitat_prop[check_indices_habprop[i],]));
  }
  for(k in 1:check_N){
    detected_prob[k] = inv_logit(detectability_bias + duration_coef*duration_minutes[check_indices_others[k]] + 
                                distance_coef*distance_traveled[check_indices_others[k]] + 
                                forest_prop_coef*forest_prop[check_indices_others[k]] + 
                                a1*cos(2*pi()*start_time[check_indices_others[k]]/24) + b1*sin(2*pi()*start_time[check_indices_others[k]]/24) + 
                                a2*cos(2*pi()*2*start_time[check_indices_others[k]]/24) + b2*sin(2*pi()*2*start_time[check_indices_others[k]]/24));
  }
  
  for(k in 1:check_N){
    observed_check[k] = bernoulli_rng(occupied_prob[k] * detected_prob[k]);
  }
}
