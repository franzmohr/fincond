#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' KPS EWMA algorithm
//' 
//' Implementation of something
//' 
//' @param y variables of the measurement equation
//' @param x variable of the factor equation
//' @param p integer of the lag order of the VAR model
//' @param const logical specifying whether a constant should be included in the VAR model
//' @param Q0 matrix of the initital values
//' @param beta0 vector of initial values
//' @param R0 matrix of the initital values
//' @param V0 matrix of the initital values
//' @param lambda0 matrix of initial values
//' @param W0 matrix of the initital values
//' @param kappaQ forgetting factor
//' @param kappaR decay factor
//' @param kappaV forgetting factor
//' @param kappaW decay factor
//' 
//' @details Produces a draw
//'
//' @noRd
//' 
// [[Rcpp::export]]
Rcpp::List kfs_ewma_parameters(arma::mat y, arma::mat x, int p, bool constant,
                               arma::mat Q0, arma::mat beta0, arma::mat R0,
                               arma::mat V0, arma::mat lambda0, arma::mat W0,
                               double kappaQ, double kappaR, double kappaV, double kappaW) {
  
  // Get dimensions
  int n_const = 0;
  if (constant) {
    n_const = 1;
  }
  int tt = y.n_cols;
  int n_y = y.n_rows;
  int n_x = x.n_rows;
  int n_beta = n_y * (n_y * p + n_const);
  //int n_lambda = n_x * n_y;
  
  // Data matrices
  arma::mat zz = arma::zeros<arma::mat>(n_y * p + n_const, tt);
  for (int i = 0; i < p; i++){
    zz.submat(i * n_y, p , (i + 1) * n_y - 1, tt - 1) = y.cols(p - i - 1, tt - 2 - i);
  }
  if (constant) {
    zz.row(n_y * p).fill(1);
  }
  arma::mat Z = arma::kron(arma::trans(zz), arma::eye<arma::mat>(n_y, n_y));
  arma::mat Z_t;
  
  arma::mat Z_lambda = arma::kron(arma::trans(y), arma::eye<arma::mat>(n_x, n_x));
  
  // Initialise matrices
  arma::mat ee = arma::zeros<arma::mat>(n_y * tt, n_y);
  arma::mat Q = arma::zeros<arma::mat>(n_y * tt, n_y);
  arma::mat beta = arma::zeros<arma::mat>(n_beta, tt);
  arma::mat beta_pred = beta0;
  arma::mat R = arma::zeros<arma::mat>(n_beta * tt, n_beta);
  arma::mat R_pred = R0;
  
  arma::mat uu = arma::zeros<arma::mat>(n_x * tt, n_x);
  arma::mat V = arma::zeros<arma::mat>(n_x * tt, n_x);
  arma::mat lambda = arma::zeros<arma::mat>(n_x * tt, n_y);
  arma::mat lambda_pred = lambda0;
  arma::cube W = arma::zeros<arma::cube>(n_y * tt, n_y, n_x);
  arma::cube W_pred = arma::zeros<arma::cube>(n_y, n_y, n_x);
  
  arma::mat F_i;
  arma::mat K;
  
  arma::vec e = arma::zeros<arma::vec>(n_y);
  ee.submat(0, 0, n_y - 1, n_y - 1) = Q0;
  Q.submat(0, 0, n_y - 1, n_y - 1) = Q0;
  beta.col(0) = beta0;
  R.submat(0, 0, n_beta - 1, n_beta - 1) = R0;
  
  arma::vec u = x.col(0) - lambda0 * y.col(0);
  uu.submat(0, 0, n_x - 1, n_x - 1) = u * arma::trans(u);
  V.submat(0, 0, n_x - 1, n_x - 1) = V0;
  lambda.rows(0, n_x - 1) = lambda0;
  for (int j = 0; j < n_x; j ++) {
    W_pred.slice(j) = W0;
    W.slice(j).rows(0, n_y - 1) = W0;
  }
  
  
  for (int j = 0; j < n_x; j++) {
    F_i = 1 / (V(j, j) + arma::trans(y.col(0)) * W.slice(j).rows(0, n_y - 1) * y.col(0));
    K = W.slice(j).rows(0, n_y - 1) * y.col(0) * F_i;
    lambda_pred.row(j) = lambda_pred.row(j) + arma::trans(K * (x(j, 0) - lambda_pred.row(j) * y.col(0)));
    W_pred.slice(j) = 1 / kappaW * (W_pred.slice(j) - K * arma::trans(y.col(0)) * W_pred.slice(j));
  } 
  
  // Kalman filter
  for (int i = 1; i < tt; i++){
    // Measurement equation
    Z_t = Z.rows(n_y * i, n_y * (i + 1) - 1);
    
    beta.col(i) = beta_pred;
    R.rows(n_beta * i, n_beta * (i + 1) - 1) = R_pred;
    
    if (i <= p) {
      ee.rows(n_y * i, n_y * (i + 1) - 1) = 0.1 * y.col(i) * arma::trans(y.col(i));
    } else {
      e = y.col(i) - Z_t * beta_pred;
      ee.rows(n_y * i, n_y * (i + 1) - 1) = e * arma::trans(e);
    }
    Q.rows(n_y * i, n_y * (i + 1) - 1) = kappaQ * Q.rows(n_y * (i - 1), n_y * i - 1) + (1 - kappaQ) * ee.rows(n_y * i, n_y * (i + 1) - 1);
    
    if (i > p) {
      F_i = arma::inv(Z_t * R_pred * arma::trans(Z_t) + Q.rows(n_y * i, n_y * (i + 1) - 1));
      K = R_pred * arma::trans(Z_t) * F_i;
      beta_pred = beta_pred + K * e;
      R_pred = 1 / kappaR * (R_pred - K * Z_t * R_pred);
    }
  
    // Factor equation
    lambda.rows(n_x * i, n_x * (i + 1) - 1) = lambda_pred;
    for (int j = 0; j < n_x; j++) {
      W.slice(j).rows(0, n_y - 1) = W_pred.slice(j);  
    }
    
    u = x.col(i) - lambda_pred * y.col(i);
    uu.rows(n_x * i, n_x * (i + 1) - 1) = u * arma::trans(u);
    V.rows(n_x * i, n_x * (i + 1) - 1) = kappaV * V.rows(n_x * (i - 1), n_x * i - 1) + (1 - kappaV) * uu.rows(n_x * i, n_x * (i + 1) - 1);
    
    for (int j = 0; j < n_x; j++) {
      F_i = 1 / (V(n_x * i + j, j) + arma::trans(y.col(i)) * W.slice(j).rows(n_y * i, n_y * (i + 1) - 1) * y.col(i));
      K = W.slice(j).rows(n_y * i, n_y * (i + 1) - 1) * y.col(i) * F_i;
      lambda_pred.row(j) = lambda_pred.row(j) + arma::trans(K * (x(j, i) - lambda_pred.row(j) * y.col(i)));
      W_pred.slice(j) = 1 / kappaW * (W_pred.slice(j) - K * arma::trans(y.col(i)) * W_pred.slice(j));
    } 
  }
  
  // Kalman smoother
  
  
  return Rcpp::List::create(Rcpp::Named("z") = lambda);
}

/*** R
load("~/fincond/data/koop.rda")
source("~/fincond/R/normalise.R")

y <- t(koop$macro)
x <- normalise(koop$financial)
x[is.na(x)] <- 0
x <- t(x)

p <- 4
tt <- ncol(y)
n_y <- nrow(y)
n_x <- nrow(x)
n_beta <- n_y * (n_y * p + 1)
n_lambda <- n_x * n_y

Q0 <- diag(1, n_y)
beta0 <- matrix(0, n_beta)
R0 <- diag(4, n_beta)

V0 <- diag(1, n_x)
lambda0 <- matrix(0, n_x, n_y)
W0 <- diag(4, n_y)

temp <- kfs_ewma_parameters(y = y, x = x, p = p, constant = TRUE,
                    Q0 = Q0, beta0 = beta0, R0 = R0,
                    V0 = V0, lambda0 = lambda0, W0 = W0,
                    kappaQ = .96, kappaR = .99, kappaV = .96, kappaW = .99)$z

temp
*/