
// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
#include <RcppEigen.h>
#include <Rcpp.h>
#include <numeric>
#include "VarParam.h"

using namespace std;
using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;

SEXP extractPriors(VarParam theta_cpp) {

	List prior;

	prior["m_0"]            = theta_cpp.m_0;
	prior["Lambda_0"]       = theta_cpp.Lambda_0;
	prior["Lambda0_m0"]     = theta_cpp.Lambda0_m0;
	prior["m0_Lambda0_m0"]  = theta_cpp.m0_Lambda0_m0;
	

	prior["a_0"]            = theta_cpp.a_0;
	prior["b_0"]            = theta_cpp.b_0;

	prior["g_0"]            = theta_cpp.g_0;
	prior["Sigma_0"]        = theta_cpp.Sigma_0;


	return wrap(prior);

	// return 0;
}



SEXP extractVarParam(VarParam theta_cpp) {

	List theta;

	theta["L"]          = theta_cpp.L;
	theta["beta_k"]     = theta_cpp.beta_k;
	theta["tau_k"]      = theta_cpp.tau_k;
	theta["gamma_k"]    = theta_cpp.gamma_k;
	theta["curr"]       = theta_cpp.curr;

	/* other variational parameters not needed to calculate approximate density,
	   but we return a few extra for debugging purposes (add as needed) */
	theta["r_nk"]       = theta_cpp.r_nk;
	
	theta["a_k"]        = theta_cpp.a_k;
	theta["b_k"]        = theta_cpp.b_k;
	
	theta["alpha"]      = theta_cpp.alpha;
	theta["xi"]         = theta_cpp.xi;
	theta["lambda"]     = theta_cpp.lambda;
	theta["phi"]        = theta_cpp.phi;

	theta["zeta_k"]     = theta_cpp.zeta_k;
	theta["m_k"]        = theta_cpp.m_k;

	theta["eta_k"]      = theta_cpp.eta_k;
	theta["mu_k"]       = theta_cpp.mu_k;

	theta["Q_k"]        = theta_cpp.Q_k;
	theta["Q_k_inv"]    = theta_cpp.Q_k_inv;

	theta["V_k"]        = theta_cpp.V_k;
	theta["V_k_inv"]    = theta_cpp.V_k_inv;

	return wrap(theta);

	// return 0;
}

// [[Rcpp::export]]
SEXP testConstructor(MAP_VEC y, MAP_MAT X, int N, int D, int K, 
			         bool intercept, int max_iter) {

	VarParam theta_cpp(y, X, N, D, K, intercept, max_iter);

	
	return extractPriors(theta_cpp);


} // end mat_list_ops() function



