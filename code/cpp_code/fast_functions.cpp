
// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>
#include <numeric>

using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using arma::vec;


#define VEC_TYPE  VectorXd
#define MAT_TYPE  MatrixXd

// #define Eigen::Map<Eigen::MatrixXd> MapMatd
// #define Eigen::Map<Eigen::VectorXd> MapVecd


// [[Rcpp::export]]
SEXP matmultC(const Eigen::Map<MAT_TYPE> A,
			  const Eigen::Map<MAT_TYPE> B) {

	MAT_TYPE C = A * B;

	return wrap(C);

} // end of matmultC() function


// [[Rcpp::export]]
SEXP quadMult(Eigen::Map<MAT_TYPE> A, Eigen::Map<VEC_TYPE> x) {

	double b = x.transpose() * (A * x);
	return wrap(b);

} // end of quadMult() function


// [[Rcpp::export]]
SEXP schur(Eigen::Map<VEC_TYPE> x, Eigen::Map<VEC_TYPE> y) {
    
    VectorXd z = x.cwiseProduct(y);
    return wrap(z);

} // end of schur() function


// [[Rcpp::export]]
SEXP testFuncs(Eigen::Map<MAT_TYPE> X_mu,      
	           Eigen::Map<VEC_TYPE> alpha, 
	           Eigen::Map<VEC_TYPE> xQx) {
    
    int n;
    int N = X_mu.rows();
    int K = xQx.size();
	MAT_TYPE res_mat(N, K); // (N x K)
	res_mat.setZero(N, K); 
	
	// fill res_mat row by row
	for (n = 0; n < N; n++) {

		// R code: sqrt((X_mu[n,] - theta$alpha[n])^2 + xQx)
		VEC_TYPE alpha_vec(K);        // K-dim vector
		alpha_vec.fill(alpha(n));
		VEC_TYPE xmu_n = X_mu.row(n);

		// Rcpp::Rcout << xmu_n.size() << std::endl;
		// Rcpp::Rcout << alpha_vec.size() << std::endl;

		VEC_TYPE vec_n = ((xmu_n - alpha_vec).array().square()).cwiseSqrt();

		res_mat.row(n) = vec_n;
	}

	return wrap(res_mat);

} // end of testFuncs() function

  
// [[Rcpp::export]]
MAT_TYPE lambda_xi_cpp(Eigen::Map<MAT_TYPE> A) {

	// note the return type is no longer SEXP since we aren't returning to R
	// 1 / (4 * xi) * tanh(0.5 * xi)
    MAT_TYPE res = (1 / (4 * A.array())).cwiseProduct((0.5 * A.array()).tanh());	

	return res;
} // end lambda_xi() function



/** Note: the input variables will be modified since they are 'Map' objects;
  *       we are writing over the objects directly
  */
// [[Rcpp::export]]
SEXP mainFunc(Eigen::Map<MAT_TYPE> lambda,       // (N x K)
	          Eigen::Map<MAT_TYPE> X,            // (N x D)
	          Eigen::Map<MAT_TYPE> X_mu,         // (N x K)
		      Eigen::Map<MAT_TYPE> Qk,           // (D x D)
		      Eigen::Map<MAT_TYPE> xi) {         

	int n, k;

	int N = lambda.rows();
	int K = lambda.cols();

	// these quantities are re-computed during every iteration of the m-step
	VEC_TYPE alpha(N);
	VEC_TYPE xQx(K);
	VEC_TYPE phi(N);


	VEC_TYPE ONES_K = VEC_TYPE::Ones(K); // K-dim vector of ones

	/** Quantities to update ---------------------------------------------------
	  * (0.1) alpha     (N x 1) -- done
	  * (0.2) xi        (N x K) -- done
	  * (0.3) phi       (N x 1) -- done
	  * (0.4) lambda    (N x K) -- done
	  ---------------------------------------------------------------------- **/
	for (n = 0; n < N; n++) {
		
		// row needs to first be stored as vec in order to do mat mults
		VEC_TYPE x_n      = X.row(n);
		VEC_TYPE xmu_n    = X_mu.row(n);         // K-dim vector
		VEC_TYPE lambda_n = lambda.row(n); 

		/* (0.1) update alpha ----------------------------------------------- */
		alpha(n) = 1 / lambda_n.sum() * 
				 (0.5 * (0.5 * K - 1) + (xmu_n.transpose() * lambda_n).value());

		for (k = 0; k < K; k++) {
			// final version of this function should have Qk changing with k
			// in the innner for loop
			xQx(k) = lambda.col(k).sum() * 
								(x_n.transpose() * (Qk * x_n)).value();

		} // end of inner for (k)
		
		/* (0.2) update xi -------------------------------------------------- */
		VEC_TYPE alpha_n(K);                     // K-dim vector
		alpha_n.fill(alpha(n));
		VEC_TYPE xi_n = ((xmu_n - alpha_n).array().square() + 
											xQx.array()).cwiseSqrt();
		xi.row(n) = xi_n;

		/* (0.3) update phi ------------------------------------------------- */
		
		VEC_TYPE log_term = (ONES_K.array() + xi_n.array().exp()).log();

		phi(n) = ((xmu_n - alpha(n) * ONES_K - xi_n) / 2 + log_term).sum(); 

	} // end of outer for (n)


	/* (0.4) update lambda (can be done for entire matrix) -------------------*/
	// Eigen::Map<MAT_TYPE> xi_map = xi;
	lambda = lambda_xi_cpp(xi);


	
	List ret;
	ret["alpha"]  = alpha;
	ret["xQx"]    = xQx;
	ret["xi"]     = xi;
	ret["phi"]    = phi;
	ret["lambda"] = lambda;

	return wrap(ret);

} // end of mainFunc()


