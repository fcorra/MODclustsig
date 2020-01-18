#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double HellDist(arma::mat X){
  double out = arma::accu(pow(X.row(0).t() - X.row(1).t(), 2));
  out = sqrt(out) / arma::datum::sqrt2;
  return out;
}
