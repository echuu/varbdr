
// VarParam.cpp -- cpp file for variational parameter object

// compile this file: g++ -I /usr/include/eigen3/ VarParam.cpp -o vb

#include "VarParam.h"
#include <stdio.h>

// VarParam constructor -- see VarParam.h for variable descriptions
VarParam::VarParam (MAP_VEC y, MAP_MAT X, int N, int D, int K, 
			        bool intercept, int max_iter) {
	

	MAT_TYPE I_D   = MAT_TYPE::Identity(D, D);

	this->L        = VEC_TYPE::Zero(max_iter);
	
	this->beta_k   = MAT_TYPE::Zero(D, K);       // (D x K) 
	this->gamma_k  = MAT_TYPE::Zero(D, K);       // (D x K)
	this->tau_k    = VEC_TYPE::Ones(K);          // (K x 1)

	this->r_nk     = MAT_TYPE::Zero(N, K);       // (N x K)
	this->log_r_nk = MAT_TYPE::Zero(N, K);       // (N x K)
	this->N_k      = VEC_TYPE::Zero(K);          // (K x 1)

	
	/* (2) var parameters for tau_k ~ Ga ( a_k, b_k )                      */
	this->a_k      = VEC_TYPE::Ones(K);          // (K x 1)
	this->b_k      = VEC_TYPE::Ones(K);          // (K x 1)


	/* (3) var parameters for gamma_k (via Bouchard)
           xi     : initialization cannot be 0 because lambda involves 1/xi
           lambda : compute here since alpha is a function of lambda	   */

	this->alpha    = VEC_TYPE::Ones(N);          // (N x 1)
	this->xi       = MAT_TYPE::Ones(N, K);       // (N x K)
	this->lambda   = lambda_xi(this->xi);        // (N x K) -- need to calculate
	this->phi      = VEC_TYPE::Zero(N);          // (N x 1)


	/* random initialization for m_k, mu_k ------------------------------- */
	// zeta_k, m_k for q(beta_k | tau_k)
	this->zeta_k   = MAT_TYPE::Zero(D, K);       // (D x K)
	this->m_k      = MAT_TYPE::Random(D, K);     // (D x K)

	// eta_k, mu_k for q(gamma_k)
	this->eta_k    = MAT_TYPE::Zero(D, K);       // (D x K)
	this->mu_k     = MAT_TYPE::Random(D, K);     // (D x K)


	/* list initialization for V_k, V_k_inv, Q_k, Q_k_inv ---------------- */
	
	/* (1) var parameters for beta_k | tau_k ~ N ( m_k, (tau_k V_k)^{-1} ) */
	/* (4) var parameters for gamma_k ~ N ( gamma_k | mu_k, Q_k^{-1} )     */
	for (int k = 0; k < K; k++) {
		this->V_k.push_back(I_D);
		this->V_k_inv.push_back(I_D);
		this->Q_k.push_back(I_D);
		this->Q_k_inv.push_back(I_D);
	}

	// list iterators
	this->Qk_it     = Q_k.begin();
	this->Qk_inv_it = Q_k_inv.begin();
	this->Vk_it     = V_k.begin();
	this->Vk_inv_it = V_k_inv.begin();


	this->curr      = 0; // current iteration of cavi

	// private data types ----------------------------------------------------

	this->N         = N;
	this->D         = D;
	this->K         = K;
	this->max_iter  = max_iter;
	this->intercept = intercept;

} // end of VarParam constructor

/* basic getters ---------------------------------------------------------- */

int VarParam::getN() {
	return this->N;
}

int VarParam::getCovDim() {
	return this->D;
}

int VarParam::getClusters() {
	return this->K;
}

int VarParam::getMaxIter() {
	return this->max_iter;
}

bool VarParam::getIntercept() {
	return this->intercept;
}


MAT_TYPE VarParam::lambda_xi (MAT_TYPE A) {
	// 1 / (4 * xi) * tanh(0.5 * xi)
    MAT_TYPE res = (1 / (4 * A.array())).cwiseProduct((0.5 * A.array()).tanh());

	return res;
} // end of lambda_xi() function


// end of VarParam.cpp file


int main(void) {
	printf("inside variational parameter class file\n");
	return 0;
}

