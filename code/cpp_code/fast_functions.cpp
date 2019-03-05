
// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>
#include <numeric>

using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using arma::vec;

// [[Rcpp::export]]
SEXP matmultC(const Eigen::Map<MatrixXd> A,
			  const Eigen::Map<MatrixXd> B) {

	MatrixXd C = A * B;

	return wrap(C);

} // end of matmultC() function


// [[Rcpp::export]]
SEXP quadMult(Eigen::Map<MatrixXd> A, Eigen::Map<VectorXd> x) {

	double b = x.transpose() * (A * x);
	return wrap(b);

} // end of quadMult() function


// [[Rcpp::export]]
SEXP schur(Eigen::Map<VectorXd> x, Eigen::Map<VectorXd> y) {
    
    VectorXd z = x.cwiseProduct(y);
    return wrap(z);

} // end of schur() function


// [[Rcpp::export]]
SEXP multC(Eigen::Map<MatrixXd> A) {
    
    VectorXd z = A.row(0);

    // VectorXd prod = A * z;

    double res = (z.transpose() * A * z);

    return wrap(res);

} // end of schur() function


// [[Rcpp::export]]
SEXP mainFunc(Eigen::Map<MatrixXd> lambda, Eigen::Map<MatrixXd> X,
		      Eigen::Map<MatrixXd> Qk) {

	int n, k;

	// rowsum of matrix
	int N = lambda.rows();
	int K = lambda.cols();

	VectorXd l_sum(N);
	VectorXd xQx(K);

	for (n = 0; n < N; n++) {

		VectorXd x_n = X.row(n);

		l_sum(n) = lambda.row(n).sum() * (x_n.transpose() * x_n).value();

		for (k = 0; k < K; k++) {

			VectorXd x_n = X.row(n);

			// xQx(k) = (X.row(n).transpose() * Qk * X.row(n)).value();
			xQx(k) = lambda.col(k).sum() * (x_n.transpose() * (Qk * x_n)).value();

		 } // end of inner for
	
	} // end of outer for

	List ret;
	ret["l_sum"] = l_sum;
	ret["xQx"]   = xQx;

	return wrap(ret);
}


