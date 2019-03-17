
// VarParam.cpp -- cpp file for variational parameter object

// compile this file: g++ -I /usr/include/eigen3/ VarParam.cpp -o vb

#include "VarParam.h"
#include <stdio.h>
#include <stdlib.h>     /* srand, rand */


// VarParam constructor -- see VarParam.h for variable descriptions
VarParam::VarParam (MAP_MAT y, MAP_MAT X, int N, int D, int K, 
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
	this->lambda   = lambda_xi(this->xi);        // (N x K)
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


	// initialize the prior parameters ---------------------------------------
	
	this->m_0           = X.colwise().mean(); // column means of X
	this->Lambda_0      = I_D;
	
	
	this->Lambda0_m0    = (this->Lambda_0) * (this->m_0);
	
	this->m0_Lambda0_m0 = ((this->m_0).transpose() * this->Lambda0_m0).value();
	

	this->a_0           = VEC_TYPE::Ones(K);
	this->b_0           = VEC_TYPE::Ones(K);

	this->g_0           = VEC_TYPE::Zero(D);
	this->Sigma_0       = I_D;

	this->y             = y;
	this->X             = X;

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
    MAT_TYPE l = (1 / (4 * A.array())).cwiseProduct((0.5 * A.array()).tanh());

	return l;
} // end of lambda_xi() function



void VarParam::mStep() {

	int n, k;

	int N = getN(); 
	int K = getClusters();

	VEC_TYPE xQx(K);                             // qty computed per iter
	VEC_TYPE ONES_K = VEC_TYPE::Ones(K);         // K-dim vector of ones

	MAT_TYPE X_mu   = (this->X) * (this->mu_k);  // (N x K) 

	/** Quantities to update ---------------------------------------------------
	  * (0.1) alpha     (N x 1) -- TODO
	  * (0.2) xi        (N x K) -- TODO
	  * (0.3) phi       (N x 1) -- TODO
	  * (0.4) lambda    (N x K) -- TODO
	  -------------------------------------------------------------------- **/

	for (n = 0; n < N; n++) {
		
		// row needs to first be stored as vec in order to do mat mults
		VEC_TYPE x_n      = this->X.row(n);      // (D x 1) : n-th row of X
		VEC_TYPE xmu_n    = X_mu.row(n);         // (K x 1) : n-th row of X_mu
		VEC_TYPE lambda_n = this->lambda.row(n); // (K x 1)

		/* (0.1) update alpha --------------------------------------------- */
		this->alpha(n) = 1 / lambda_n.sum() * 
			   (0.5 * (0.5 * K - 1) + (xmu_n.transpose() * lambda_n).value());

		for (k = 0; k < K; k++) {
			/** final version of this function should have Qk changing with k
			    in the innner for loop **/

			xQx(k) = (x_n.transpose() * ((*this->Qk_inv_it) * x_n)).value();

			advance(this->Qk_inv_it, 1);

		} // end of inner for (k)
		
		this->Qk_inv_it = Q_k_inv.begin(); // bring back to front of the list

		/* (0.2) update xi ------------------------------------------------ */
		VEC_TYPE alpha_n(K);                     // K-dim vector
		alpha_n.fill(this->alpha(n));
		VEC_TYPE xi_n = ((xmu_n - alpha_n).array().square() + 
											xQx.array()).cwiseSqrt();
		this->xi.row(n) = xi_n;

		/* (0.3) update phi ----------------------------------------------- */
		
		VEC_TYPE log_term = (ONES_K.array() + xi_n.array().exp()).log();

		this->phi(n) = ((xmu_n - alpha(n) * ONES_K - xi_n) / 2 + log_term).sum(); 

	} // end of outer for (n)


	/* (0.4) update lambda (can be done for entire matrix) -----------------*/
	// Eigen::Map<MAT_TYPE> xi_map = xi;
	this->lambda = lambda_xi(this->xi);


} // end mStep() function



// end of VarParam.cpp file

int main(void) {
	printf("inside variational parameter class file\n");
	return 0;
}

