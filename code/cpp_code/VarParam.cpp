
// VarParam.cpp -- cpp file for variational parameter object


#include VarParam.h


// VarParam constructor -- see VarParam.h for variable descriptions
VarParam::VarParam (MAP_VEC y, MAP_MAT X, int N, int D, int K, 
			        bool intercept, int max_iter) {
	

	MAT_TYPE   I_D = MAT_TYPE::Identity(D, D);
	VEC_TYPE   L (max_iter);
	
	MAT_TYPE   beta_k   = MAT_TYPE::Zero(D, K);  // (D x K) 
	MAT_TYPE   gamma_k  = MAT_TYPE::Zero(D, K);  // (D x K)
	MAT_TYPE   r_nk     = MAT_TYPE::Zero(N, K);  // (N x K)

	MAT_TYPE   r_nk     = MAT_TYPE::Zero(N, K);  // (N x K)
	MAT_TYPE   r_nk     = MAT_TYPE::Zero(N, K);  // (N x K)
	VEC_TYPE   N_k      = VEC_TYPE::Zero(K);     // (K x 1)

	/* (1) var parameters for beta_k | tau_k ~ N ( m_k, (tau_k V_k)^{-1} ) */
	
	
	/* (2) var parameters for tau_k ~ Ga ( a_k, b_k )                      */
	

	/* (3) var parameters for gamma_k (via Bouchard)
           xi     : initialization cannot be 0 because lambda involves 1/xi
           lambda : compute here since alpha is a function of lambda	   */


	/* (4) var parameters for gamma_k ~ N ( gamma_k | mu_k, Q_k^{-1} )     */



	/* random initialization for m_k, mu_k ------------------------------- */



	// private data types ----------------------------------------------------

	this.N        = N;
	this.D        = D;
	this.K        = K;
	this.max_iter = max_iter;
	this.intercept = intercept;

} // end of VarParam constructor





/* basic getters ---------------------------------------------------------- */

int VarParam::getN() {
	return this.N;
}

int VarParam::getCovDim() {
	return this.D;
}

int VarParam::getClusters() {
	return this.K;
}

int VarParam::getMaxIter() {
	return this.max_iter;
}

bool VarParam::getIntercept() {
	return this.intercept;
}

// end of VarParam.cpp file
