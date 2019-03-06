
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
SEXP testFuncs(Eigen::Map<MatrixXd> X_mu,      
	           Eigen::Map<VectorXd> alpha, 
	           Eigen::Map<VectorXd> xQx) {
    
    int n;
    int N = X_mu.rows();
    int K = xQx.size();
	MatrixXd res_mat(N, K); // (N x K)
	res_mat.setZero(N, K); 
	
	// fill res_mat row by row
	for (n = 0; n < N; n++) {

		// R code: sqrt((X_mu[n,] - theta$alpha[n])^2 + xQx)
		VectorXd alpha_vec(K);        // K-dim vector
		alpha_vec.fill(alpha(n));
		VectorXd xmu_n = X_mu.row(n);

		// Rcpp::Rcout << xmu_n.size() << std::endl;
		// Rcpp::Rcout << alpha_vec.size() << std::endl;

		VectorXd vec_n = ((xmu_n - alpha_vec).array().square()).cwiseSqrt();

		res_mat.row(n) = vec_n;
	}

	return wrap(res_mat);

} // end of schur() function


// [[Rcpp::export]]
SEXP lambda_func(Eigen::Map<MatrixXd> A, int n) {

	VectorXd x(A.cols());      // K-dim vector
	x.fill(0);

	A.row(n) = x;

	return wrap(A);
}



// [[Rcpp::export]]
SEXP mainFunc(Eigen::Map<MatrixXd> lambda,       // (N x K)
	          Eigen::Map<MatrixXd> X,            // (N x D)
	          Eigen::Map<MatrixXd> X_mu,         // (N x K)
		      Eigen::Map<MatrixXd> Qk) {         // (D x D)

	int n, k;

	int N = lambda.rows();
	int K = lambda.cols();

	VectorXd alpha(N);
	VectorXd xQx(K);
	MatrixXd xi(N, K);


	/** Quantities to update ---------------------------------------------------
	  * (0.1) alpha     (N x 1) -- done
	  * (0.2) xi        (N x K) -- done
	  * (0.3) lambda    (N x K) -- done
	  * (0.4) phi       (N x 1) -- TODO
	  ---------------------------------------------------------------------- **/
	for (n = 0; n < N; n++) {
		// row needs to first be stored as vec
		VectorXd x_n      = X.row(n);
		VectorXd xmu_n    = X_mu.row(n);         // K-dim vector
		VectorXd lambda_n = lambda.row(n); 

		// (0.1) update alpha
		alpha(n) = lambda_n.sum() * (xmu_n.transpose() * lambda_n).value();

		for (k = 0; k < K; k++) {
			xQx(k) = lambda.col(k).sum() * 
								(x_n.transpose() * (Qk * x_n)).value();

		} // end of inner for
		
		// (0.2) update xi
		VectorXd alpha_n(K);                     // K-dim vector
		alpha_n.fill(alpha(n));
		VectorXd vec_n = ((xmu_n - alpha_n).array().square() + 
											xQx.array()).cwiseSqrt();
		xi.row(n) = vec_n;

	} // end of outer for


	// (0.3) update lambda (can be done for entire matrix)

	List ret;
	ret["alpha"] = alpha;
	ret["xQx"]   = xQx;
	ret["xi"]    = xi;

	return wrap(ret);

} // end of mainFunc()


