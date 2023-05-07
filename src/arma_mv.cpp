// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//' Matrix Product in RcppArmadillo.
//'
//' @param m numeric matrix
//' @param v numeric vector
//' @return matrix product of m and v
//' @export
// [[Rcpp::export(arma_mv)]]
arma::mat arma_mv(const arma::mat& m, const arma::vec& v) {
  return m * v;
};