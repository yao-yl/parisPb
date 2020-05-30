 functions {
   matrix cov_period(real[] x, real rho) {
     int N = size(x);
     matrix[N, N] K;
     for (i in 1:N)
        for (j in 1:N){
          K[i,j] =  exp(-2*(sin( fabs(x[i] - x[j]) / 2)^2) / (rho ^ 2));
        }
     return K;
     }
    matrix cov_period_across(real[] x1, real[] x2, real rho) {
     int N1 = size(x1);
     int N2 = size(x2);
     matrix[N1, N2] K;
     for (i in 1:N1)
        for (j in 1:N2){
          K[i,j] =  exp(-2*(sin( fabs(x1[i] - x2[j]) / 2)^2) / (rho ^ 2));
        }
     return K;
     } 
                       
  vector gp_pred_rng(real[] d2, real[] theta2,
                     vector y1, real[] d1, real[] theta1,
                     real alpha, real rho_d, real sigma, real rho_theta, real delta) {
    int N1 = rows(y1);
    int N2 = size(d2);
    vector[N2] f2;
    {
      matrix[N1, N1] K =   cov_exp_quad(d1, alpha, rho_d)  .* cov_period(theta1, rho_theta)
                           + diag_matrix(rep_vector(square(sigma), N1)); 
      matrix[N1, N1] L_K = cholesky_decompose(K);
      vector[N1] L_K_div_y1 = mdivide_left_tri_low(L_K, y1);
      vector[N1] K_div_y1 = mdivide_right_tri_low(L_K_div_y1', L_K)';
      matrix[N1, N2] k_x1_x2 = cov_exp_quad(d1, d2, alpha, rho_d)  .* cov_period_across(theta1,theta2,rho_theta);
      vector[N2] f2_mu = (k_x1_x2' * K_div_y1);
      matrix[N1, N2] v_pred = mdivide_left_tri_low(L_K, k_x1_x2);
      matrix[N2, N2] cov_f2 =   cov_exp_quad(d2, alpha, rho_d) .* cov_period(theta2, rho_theta)  - v_pred' * v_pred
                              + diag_matrix(rep_vector(delta, N2));
      f2 = multi_normal_rng(f2_mu, cov_f2);
    }
    return f2;
  }
}

data {
  int<lower=1> N;
  real theta[N];
  real d[N];
  vector[N] y;
  int type[N];
  int<lower=1> N_test;
  real d_test[N_test];
  real theta_test[N_test];
}

parameters {
  real<lower=0> rho_d;
  real<lower=0> rho_theta;
  real<lower=0> alpha;
  real<lower=0> sigma;
  vector[4] type_eff_short; 
}
transformed parameters {
  vector[5] type_eff; 
  type_eff[1:4]=type_eff_short;
  type_eff[5]=-sum(type_eff_short);
}
model {
  matrix[N, N] cov =   cov_exp_quad(d, alpha, rho_d) .* cov_period(theta, rho_theta)
  + diag_matrix(rep_vector(square(sigma), N));
  matrix[N, N] L_cov = cholesky_decompose(cov);
  y ~ multi_normal_cholesky(type_eff[type] , L_cov);
  rho_d ~ normal(0,1.5);
  rho_theta ~ normal(0,1);
  alpha ~ normal(0,6);
  sigma~ normal(0,6);
  type_eff_short ~ normal(0,1);
}

generated quantities {
  vector[N_test] f_predict = gp_pred_rng(d_test, theta_test, y,  d, theta,
                      alpha,  rho_d,  sigma,  rho_theta,  1e-9);
   vector[N_test] y_predict;
   {
       vector[N_test] noise= to_vector (normal_rng ( rep_vector(0,N_test),  rep_vector(1,N_test) ));
       y_predict=f_predict + noise* sigma;
   }
}
