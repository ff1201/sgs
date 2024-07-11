// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//' Matrix Product in RcppArmadillo.
//'
//' @param m numeric sparse matrix
//' @param v numeric vector
//' @return matrix product of m and v
//' @export
// [[Rcpp::export(arma_sparse)]]
arma::mat arma_sparse(const arma::sp_mat& m, const arma::vec& v) {
  return m * v;
};