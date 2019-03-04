// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;


// [[Rcpp::export]]
SEXP armaMatMult(arma::mat A, arma::mat B){
    arma::mat C = A * B;

    return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP eigenMatMult(Eigen::MatrixXd A, Eigen::MatrixXd B){
    Eigen::MatrixXd C = A * B;

    return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP eigenMapMatMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
    Eigen::MatrixXd C = A * B;

    return Rcpp::wrap(C);
}


// [[Rcpp::export]]
SEXP schur(arma::vec x, arma::vec y){
    arma::vec z = x % y;

    return Rcpp::wrap(z);
}


// [[Rcpp::export]]
NumericVector fastRowSum(NumericMatrix x) {
	
	int nrow = x.nrow(), ncol = x.ncol();
	
	NumericVector out(nrow);

	for (int i = 0; i < nrow; i++) {
		double total = 0;
	    for (int j = 0; j < ncol; j++) {
	    	total += x(i, j);
	    }
	    out[i] = total;
	}
	
	return out;
}

