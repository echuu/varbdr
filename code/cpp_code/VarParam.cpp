
// VarParam.cpp -- cpp file for variational parameter object

// compile this file: g++ -I /usr/include/eigen3/ VarParam.cpp -o vb

#include "VarParam.h"
#include <stdio.h>
#include <stdlib.h>     /* srand, rand */


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
	int D = getCovDim();
	int K = getClusters();

	VEC_TYPE xQx(K);                             // K-dim vector used in xi
	VEC_TYPE ONES_K = VEC_TYPE::Ones(K);         // K-dim vector of ones
	VEC_TYPE ONES_N = VEC_TYPE::Ones(N);         // N-dim vector of ones
	MAT_TYPE I_D   = MAT_TYPE::Identity(D, D);   // (D x D) identity matrix
 
	MAT_TYPE X_mu   = (this->X) * (this->mu_k);  // (N x K) : X * mu_k

	MAT_TYPE rl_nk_xx;                           // used in q(gamma)
	MAT_TYPE r_nk_xx; 							 // used in q(beta|tau)
	MAT_TYPE r_x;                                // used in rl_nk_xx, r_nk_xx

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


	/* (0.4) update lambda (can be done for entire matrix) ----------------- */
	// Eigen::Map<MAT_TYPE> xi_map = xi;
	this->lambda = lambda_xi(this->xi);


	/** Quantities to update -------------------------------------------------
	  * (1) update variational distributions/parameters
	  *     (1.1) q(gamma)     DONE
	  *     (1.2) q(beta|tau)  DONE
	  *     (1.3) q(tau)       TODO
	  * 
	  * (2) update posterior means
	  *     (2.1) beta_k       TODO
	  *     (2.2) tau_k        TODO
	  *     (2.3) gamma_k      TODO
	  ----------------------------------------------------------------------  */

	VEC_TYPE y2  = y.array().square(); // move this into constructor later since
									   // y is fixed, save calculation

	/* q(gamma) : Q_k, Q_k^{-1}, eta_k, mu_k -------------------------------- */
	for (k = 0; k < K; k++) {

		rl_nk_xx = MAT_TYPE::Zero(D, D); // reset entries to 0 every iter
		r_nk_xx  = MAT_TYPE::Zero(D, D); // reset entries to 0 every iter

		// calculate intermediate quantities for q(gamma), q(beta|tau)
		for (n = 0; n < N; n++) {
			VEC_TYPE x_n = this->X.row(n);
			// x_n * x_n' can be precomputed 
			r_x      = this->r_nk(n,k) * (x_n * x_n.transpose()).array();
			rl_nk_xx = rl_nk_xx.array() + this->lambda(n, k) * r_x.array(); 
			r_nk_xx  = r_nk_xx + r_x;

		} // end inner for

		/* (1.1) q(gamma)    : Q_k, Q_k_inv, eta_k, mu_k -------------------- */
		
		VEC_TYPE rnk_k = this->r_nk.col(k);                        // (N x 1)
		VEC_TYPE lam_k = this->lambda.col(k);                      // (N x 1)
 		VEC_TYPE alpha_lam_k = 0.5 * ONES_N + 2 * lam_k.cwiseProduct(this->alpha);  

		*this->Qk_it      = I_D.array() + 2 * rl_nk_xx.array();
		*this->Qk_inv_it  = (*this->Qk_it).inverse();
		
		this->eta_k.col(k) = this->X.transpose() * (rnk_k.cwiseProduct(alpha_lam_k)); // (D x 1)
		this->mu_k.col(k)  = (*this->Qk_inv_it) * this->eta_k.col(k);

		/* (1.2) q(beta|tau) : V_k, V_k_inv, zeta_k, m_k -------------------- */
		VEC_TYPE rnk_y = rnk_k.cwiseProduct(this->y);

		*this->Vk_it     = this->Lambda_0 + r_nk_xx;
		*this->Vk_inv_it = (*this->Vk_it).inverse();
		
		this->zeta_k.col(k) = this->Lambda0_m0.array() + (this->X.transpose() * rnk_y).array();
		this->m_k.col(k)    = (*this->Vk_inv_it) * this->zeta_k.col(k);

		// advance precision matrix iterators
		advance(this->Qk_it, 1);
		advance(this->Qk_inv_it, 1);
		advance(this->Vk_it, 1);
		advance(this->Vk_inv_it, 1);

	} // end of outer for -- finish q(gamma), q(beta|tau) updates


	// reset iterators
	this->Qk_it     = Q_k.begin(); 
	this->Qk_inv_it = Q_k_inv.begin(); 
	this->Vk_it     = V_k.begin(); 
	this->Vk_inv_it = V_k_inv.begin(); 


	/** IN PROGRESS -------------------------------------------------------- **/

	/* (1.3) update q(tau)   : a_k, b_k ------------------------------------- */



	/* (2.1) update beta_k  ------------------------------------------------- */


	/* (2.2) update tau_k   ------------------------------------------------- */


	/* (2.3) update gamma_k ------------------------------------------------- */


	// update current iteration
	this->curr++;



} // end mStep() function



// end of VarParam.cpp file

int main(void) {
	printf("inside variational parameter class file\n");
	return 0;
}

