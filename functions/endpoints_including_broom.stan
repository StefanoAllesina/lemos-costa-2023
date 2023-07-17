functions {
  vector x_k(vector presence, matrix A, int k) {
    int n = num_elements(presence);
    int indices[k];
    vector[k] onesk;
    matrix[k,k] Ak;
    int j = 1;
    for (i in 1:n){
      if (presence[i] > 0.0) {
        indices[j] = i;
        j += 1;
      }
    }
    for (i in 1:k){
      onesk[i] = 1.0;
      for (l in 1:k){
        Ak[i,l] = A[indices[i], indices[l]];
      }
    }
    return Ak \ onesk;
  }
  matrix get_endpoints(matrix P, matrix A, int[] nsp){
    // P is transposed: 
    // each column is one experiment
    int m = cols(P);
    int n = rows(P);
    int l;
    matrix[n, m] X;
    for (i in 1:m){
      vector[nsp[i]] xtmp;
      xtmp = x_k(col(P, i), A, nsp[i]);
      l = 1;
      for (j in 1:n){
        if (P[j, i] > 0.0){
          X[j, i] = xtmp[l];
          l += 1;
        } else { 
          X[j, i] = 0.0;
        }
      }
    }
    return X;
  }
  matrix build_B(vector pars, matrix V, int modelnum){
    int n = cols(V);
    int m = rows(V);
    int np = num_elements(pars);
    vector[m] lambda;
    vector[n] theta;
    vector[n] gamma;
    matrix[m, m] L;
    matrix[n, n] Theta;
    matrix[n, n] Gamma;
    lambda = pars[1:m];
    if (modelnum == 1){
      Theta = identity_matrix(n);
      Gamma = identity_matrix(n);
    }
    if (modelnum == 2){
      gamma = fabs(pars[(m + 1):(np - 1)]);
      gamma = n * gamma / sum(gamma);
      Gamma = diag_matrix(gamma);
      Theta = identity_matrix(n);
    }
    if (modelnum == 3){
      theta = fabs(pars[(m + 1):(np - 1)]);
      theta = n * theta / sum(theta);
      Theta = diag_matrix(theta);
      Gamma = identity_matrix(n);
    }
    if (modelnum == 4){
      gamma = fabs(pars[(m + 1):(m + n)]);
      gamma = n * gamma / sum(gamma);
      Gamma = diag_matrix(gamma);
      theta = fabs(pars[(m + n + 1):(m + 2 * n)]);
      theta = n * theta / sum(theta);
      Theta = diag_matrix(theta);
    }
    L = diag_matrix(lambda);
    return Gamma * V' * L * V * Theta;
  }
  matrix get_prediction(vector pars, matrix V, matrix P, matrix Obs, int[] replicates, int[] nsp, int modelnum){
    int m = cols(P);
    int n = rows(P);
    matrix[n, m] Pred =  get_endpoints(P, build_B(pars, V, modelnum), nsp);
    int m2 = cols(Obs);
    matrix[n, m2] Pred2;
    for (i in 1:n){
      for(j in 1:m2){
        Pred2[i,j] = Pred[i, replicates[j]];
      }
    }
    return Pred2;
  }
  real likely_gamma(vector pars, matrix V, matrix P, matrix Obs, int[] replicates, int[] nsp, int modelnum){
    int m = cols(Obs);
    int n = rows(Obs);
    int np = num_elements(pars);
    matrix[n, m] X = get_prediction(pars, V, P, Obs, replicates, nsp, modelnum);
    real alpha;
    alpha = 1.0 + fabs(pars[np]);
    real likely = 0.0;
    real PENAL = -10000;
    for (i in 1:n){
      for (j in 1:m){
        if (Obs[i,j] > 0.0){
          if (X[i,j] > 0) {
            likely += gamma_lpdf(Obs[i,j] | alpha, alpha / X[i,j]);
          } else {
            likely += PENAL + gamma_lpdf(Obs[i,j] | alpha, alpha * 10000.0) - fabs(X[i,j]);
          }
        }
      }
    }
    return likely;
  }
}

