// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>

using namespace Rcpp;

//' Reimplementation of quadratic form in C++
//'
//' @param x the quadratic form
//' @param y the matrix on which x is applied
//'
//' @return t(x)%*% y %*% x, but faster I hope
// [[Rcpp::export]]
arma::mat tquadform_cpp(mat &x, mat &y) { return (x.t() * y * x); }