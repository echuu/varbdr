
// VarParam.h -- header file for variational parameter object

// #include <RcppEigen.h>
#include <numeric>
#include <Eigen/Dense>
#include <list>
//#include <Rmath.h>
#include <math.h>
#include <unsupported/Eigen/SpecialFunctions> // digamma(), lgamma()


using namespace std;
// using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;

#define VEC_TYPE  VectorXd
#define MAT_TYPE  MatrixXd
#define MAP_MAT   Eigen::Map<Eigen::MatrixXd> 
#define MAP_VEC   Eigen::Map<Eigen::VectorXd> 


/* LIST OF VARIATIONAL PARAMETERS TO INCLUDE IN THE CLASS ----------------------
 * N 	           :  # of observations
 * D    	       :  dimension of covariates
 * K    	       :  # of clusters
 * max_iter        :  maximum number of iterations to run CAVI
 * intercept	   :  boolean for inclusion of intercept (col of 1's) in X
 * L               :  vector to store the variational lower bound
 * 
 * beta_k          :  (D x K) k-th col is coeff vector for k-th gaussian 
 * tau_k           :  (K x 1) k-th elem is precision for k-th gaussian
 * gamma_k         :  (D x K) k-th col is coeff vector for k-th cluster
 *
 * r_nk            :  (N x K) matrix of responsibilities
 * log_r_nk        :  (N x K) log(r_nk)
 * N_k             :  (K x 1) column sums of r_nk
 * 
 * V_k             :  K-dim list of (D x D) inverse prec. matrices for beta_k
 * V_k_inv         :  K-dim list of (D x D) precision matrices for beta_k
 * zeta_k          :  (D x K) matrix of intermediate values 
 * m_k             :  (D x K) k-th col is mean of beta_k
 * 
 * alpha           :  (N x 1) extra variational parameter for approx E[softmax]
 * xi              :  (N x K) extra variational parameter for approx E[softmax]
 * lambda          :  (N x K) lambda(xi) -- include file w/ this func
 * phi             :  (N x 1) function of extra variational parameters
 * 
 * Q_k             :  K-dim list of (D x D) inverse prec. matrices for gamma_k
 * Q_k_inv         :  K-dim list of (D x D) precision matrices for gamma_k
 * eta_k           :  (D x K) matrix of intermediate values 
 * mu_k            :  (D x K) k-th col is mean of gamma_k
 * 
 * a_k             : (K x 1) vector of shape params for tau_k
 * b_k             : (K x 1) vector of rate params for tau_k
 ---------------------------------------------------------------------------- */


class VarParam {

	public:  // all things to be modified -> make public
		VEC_TYPE        L;
		
		MAT_TYPE        beta_k, gamma_k; 
		VEC_TYPE        tau_k;
		
		MAT_TYPE        r_nk, log_r_nk;
		VEC_TYPE        N_k;

		list<MAT_TYPE>  V_k, V_k_inv, Q_k, Q_k_inv;
		
		list<MAT_TYPE>::iterator Vk_it, Vk_inv_it, Qk_it, Qk_inv_it;

		MAT_TYPE        zeta_k, m_k, eta_k, mu_k;

		VEC_TYPE        alpha, phi;
		MAT_TYPE        xi, lambda;

		VEC_TYPE        a_k, b_k;

		int             curr;

       	/* make these private later */
       	// priors are assumed to be same for k = 1,...,K
		VEC_TYPE   m_0;           // (D x 1) prior mean for beta_k
		MAT_TYPE   Lambda_0;      // (D x D) prior covariance for beta_k
		MAT_TYPE   Lambda0_m0;    // (D x 1) Lambda_0 * m_0
		double     m0_Lambda0_m0; // (1 x 1) m_0' * Lambda_0 * m_0
		double     lgd_Lambda_0;  // (1 x 1) log(det(Lambda_0))

		VEC_TYPE   a_0, b_0;      // (K x 1) prior shape, rate params for tau_k
		VEC_TYPE   lg_a0, log_b0; // (K x 1) lgamma(a_0), log(b_0)


		MAT_TYPE   g_0;      // (D x 1) prior mean for gamma_k
		MAT_TYPE   Sigma_0;  // (D x D) prior covariance for gamma_k


		/* --- data --- */
		VEC_TYPE    y;        // (N x 1) response
		VEC_TYPE    y2;       // (N x 1) element-wise square of the response
		MAT_TYPE    X;        // (N x D) design matrix


		list<MAT_TYPE> xn_mat;   // N-dim list of (D x D) matricies 
		                                     // made of X[n,] X[n,]'

		list<MAT_TYPE>::iterator xn_mat_it;  // iterator for xn_mat


		double tol; // tolerance for evaluating convergence of elbo


		// ---------------------------------------------------------------------

		// constructor for variational parameters
		VarParam (MAP_VEC y, MAP_MAT X, int N, int D, int K, 
			      bool intercept, int max_iter);

		static MAT_TYPE lambda_xi (MAT_TYPE A);
		static VEC_TYPE my_digamma(VEC_TYPE x);

		static double   lse(VEC_TYPE vec, int dim);
		static VEC_TYPE lse_rows(MAT_TYPE X, int N, int D);

		// SEXP   extractVarparams();

		/* none of the VB routines should require any input since all 
		   calculations use internal variables                        */

		// e-step, m-step will *internally* modify variational parameters
		void       eStep();
		void       mStep();

		// elbo will *internally* compute the variational lower bound
		void     elbo();
		bool     elboConverge();

		void     nextIter();

		// basic getters
		int        getN();
		int        getCovDim();
		int        getClusters();
		int        getMaxIter();
		bool       getIntercept();

	private: // include the priors here
		int        N, D, K, max_iter;
		bool       intercept;



};
