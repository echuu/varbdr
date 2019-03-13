


// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>
#include <numeric>

using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;

#define VEC_TYPE  VectorXd
#define MAT_TYPE  MatrixXd

#define MAP_MAT Eigen::Map<Eigen::MatrixXd> 
#define MAP_VEC Eigen::Map<Eigen::VectorXd> 


// [[Rcpp::export]]
SEXP updateVariational(MAP_MAT X,                // (N x K)
					   MAP_VEC y,                // (N x 1)
	          		   MAP_MAT r_nk,             // (N x D)
	          		   MAP_VEC N_k,              // (K x 1)
	          		   MAP_MAT lambda,           // (N x K)
	          		   MAP_VEC alpha,            // (N x 1)
	          		   MAP_MAT Lambda_0,         // (D x D) 
	          		   MAP_VEC Lambda0_m0,       // (D x 1)
	          		   double  m0_Lambda0_m0,    // (1 x 1)
	          		   int     a_0,				 // (1 x 1)
	          		   int     b_0) {     		 // (1 x 1)

	int n, k;

	int N = lambda.rows();
	int K = lambda.cols();
	int D = X.cols();

	MAT_TYPE I_D = MAT_TYPE::Identity(D, D);     // (D x D) identity matrix
	VEC_TYPE y2  = y.array().square();
												 // needs to be K x (D x D)
	
	VEC_TYPE ONES_K = VEC_TYPE::Ones(K);         // K-dim vector of ones
	VEC_TYPE ONES_N = VEC_TYPE::Ones(N);         // N-dim vector of ones

	MAT_TYPE rl_nk_xx;                           // used to compute q(gamma)
	MAT_TYPE r_nk_xx; 							 // used to compute q(beta|tau)
	MAT_TYPE r_x;                                // used in rl_nk_xx, r_nk_xx


	/* -- variables below are ones that will be stored in 'theta' later ----  */
	MAT_TYPE Q_k, Q_k_inv;
	MAT_TYPE eta_k, mu_k; 

	MAT_TYPE V_k, V_k_inv;
	MAT_TYPE zeta_k, m_k; 

	VEC_TYPE a_k(K);
	VEC_TYPE b_k(K);

	MAT_TYPE beta_k, gamma_k;
	VEC_TYPE tau_k;


	/** Quantities to update ---------------------------------------------------
	  * (1) update variational distributions/parameters
	  *     (1.1) q(gamma)     IN-PROGRESS
	  *     (1.2) q(beta|tau)  IN-PROGRESS
	  *     (1.3) q(tau)       TODO
	  * 
	  * (2) update posterior means
	  *     (2.1) beta_k       TODO
	  *     (2.2) tau_k        TODO
	  *     (2.3) gamma_k      TODO
	  ----------------------------------------------------------------------  */

	/* q(gamma) : Q_k, Q_k^{-1}, eta_k, mu_k -------------------------------- */
	for (k = 0; k < K; k++) {

		rl_nk_xx = MAT_TYPE::Zero(D, D);
		r_nk_xx  = MAT_TYPE::Zero(D, D);

		// calculate intermediate quantities for q(gamma), q(beta|tau)
		for (n = 0; n < N; n++) {
			VEC_TYPE x_n = X.row(n);
			
			r_x      = r_nk(n,k) * (x_n * x_n.transpose()).array();
			rl_nk_xx = rl_nk_xx.array() + lambda(n, k) * r_x.array(); 
			r_nk_xx  = r_nk_xx + r_x;

		} // end inner for

		/* (1.1) q(gamma)    : Q_k, Q_k_inv, eta_k, mu_k -------------------- */
		
		VEC_TYPE rnk_k = r_nk.col(k);                                // (N x 1)
		VEC_TYPE lam_k = lambda.col(k);                              // (N x 1)
 		VEC_TYPE alpha_lam_k = 0.5 * ONES_N + 2 * lam_k.cwiseProduct(alpha);  

		Q_k     = I_D.array() + 2 * rl_nk_xx.array();
		Q_k_inv = Q_k.inverse();
		eta_k   = X.transpose() * (rnk_k.cwiseProduct(alpha_lam_k)); // (D x 1)
		mu_k    = Q_k_inv * eta_k;

		/* (1.2) q(beta|tau) : V_k, V_k_inv, zeta_k, m_k -------------------- */
		VEC_TYPE rnk_y = rnk_k.cwiseProduct(y);

		V_k     = Lambda_0 + r_nk_xx;
		V_k_inv = V_k.inverse();
		zeta_k  = Lambda0_m0.array() + (X.transpose() * rnk_y).array();
		m_k     = V_k_inv * zeta_k;

	} // end of outer for -- finish q(gamma), q(beta|tau) updates


	/* (1.3) update q(tau)   : a_k, b_k ------------------------------------- */
	a_k = (a_0 * ONES_K).array() + (0.5 * N_k).array();

	for (k = 0; k < K; k++) {
		VEC_TYPE rnk_k = r_nk.col(k);                                // (N x 1)
		b_k(k) = - (zeta_k.transpose() * (V_k_inv * zeta_k)).value() + 
					(rnk_k.cwiseProduct(y2)).sum();
	}

	b_k = (b_0 * ONES_K).array() + 
				(0.5 * (b_k + m0_Lambda0_m0 * ONES_K)).array();

	/* (2.1) update beta_k  ------------------------------------------------- */
	beta_k = m_k;			

	/* (2.2) update tau_k   ------------------------------------------------- */
	tau_k = a_k.cwiseQuotient(b_k);

	/* (2.3) update gamma_k ------------------------------------------------- */
	gamma_k = mu_k;

	// update current iteration

	List ret;
	ret["r_x"]         = r_x;
	ret["rl_nk_xx"]    = rl_nk_xx;
	ret["r_nk_xx"]     = r_nk_xx;

	ret["Q_k"]         = Q_k;
	ret["Q_k_inv"]     = Q_k_inv;
	ret["eta_k"]       = eta_k;
	ret["mu_k"]        = mu_k;

	ret["V_k"]         = V_k;
	ret["V_k_inv"]     = V_k_inv;
	ret["zeta_k"]      = zeta_k;
	ret["m_k"]         = m_k;

	
	ret["a_k"]         = a_k;
	ret["b_k"]         = b_k;
	ret["beta_k"]      = beta_k;
	ret["tau_k"]       = tau_k;
	ret["gamma_k"]     = gamma_k;

	return wrap(ret);

} // end updateVariational() function

