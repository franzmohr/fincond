#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' KPS EWMA algorithm
//' 
//' Implementation of something
//' 
//' @param x a time series.
//' 
//' @details Produces a draw
//'
//' @noRd
//' 
// [[Rcpp::export]]
Rcpp::List pca(arma::mat x) {
  
  arma::mat xx = arma::trans(x) * x;
  
  arma::mat coeff;
  arma::mat score;
  arma::vec latent;
  arma::vec tsquared;
  
  princomp(coeff, score, latent, tsquared, xx);
  
  return Rcpp::List::create(Rcpp::Named("coeff") = coeff,
                            Rcpp::Named("score") = score,
                            Rcpp::Named("latent") = latent,
                            Rcpp::Named("tsquared") = tsquared);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
#timesTwo(42)
*/
