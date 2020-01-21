#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat csam(arma::mat& x, arma::mat& em) {
  int ns = x.n_rows;
  int nend = em.n_rows;
  arma::mat out(ns, nend);
  out.zeros();
  arma::mat normEM = sum(pow(em,2),1);
  arma::mat normT  = sum(pow(x,2),1);
  
  for(int j = 0; j < nend; ++j){
    for(int i = 0; i < ns; ++i){
      out(i,j) = acos(accu(x.row(i) % em.row(j)) / sqrt(normT(i,0) * normEM(j,0)));
    }
  }
  return out;
}