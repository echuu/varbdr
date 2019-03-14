
// VarParam.h -- header file for variational parameter object

#include <RcppEigen.h>
#include <numeric>

using namespace std;
using namespace Rcpp;
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
		MAT_TYPE        zeta_k, m_k, eta_k, mu_k;

		VEC_TYPE        alpha, phi;
		MAT_TYPE        xi, lambda;

		VEC_TYPE        a_k, b_k;

	private: // include the priors here
		int  N, D, K, max_iter;
		bool intercept;
		// other priors

};
