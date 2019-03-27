
// VarParam.cpp -- cpp file for variational parameter object

// compile this file: g++ -I /usr/include/eigen3/ VarParam.cpp -o vb

#include "VarParam.h"
#include <stdio.h>
#include <stdlib.h>     /* srand, rand */
#include <iostream>
#include <iomanip>
#include <unsupported/Eigen/SpecialFunctions> // digamma(), lgamma()


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


	this->curr      = 1; // current iteration of cavi

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
	
	this->lgd_Lambda_0  = log(Lambda_0.determinant());

	this->a_0           = VEC_TYPE::Ones(K);
	this->b_0           = VEC_TYPE::Ones(K);

	this->lg_a0         = (this->a_0).array().lgamma();
	this->log_b0        = (this->b_0).array().log();

	this->g_0           = VEC_TYPE::Zero(D);
	this->Sigma_0       = I_D;

	this->y             = y;
	this->y2            = this->y.array().square();
	this->X             = X;



	// precompute each X[n,] X[n,]'
	for (int n = 0; n < N; n++) {
		VEC_TYPE x_n = this->X.row(n);
		this->xn_mat.push_back(x_n * x_n.transpose());
	}
	this->xn_mat_it = xn_mat.begin();

	this->tol = std::pow(10, -3);

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



double VarParam::lse(VEC_TYPE vec, int dim) {
	double maxVal = vec.maxCoeff(); 

	VEC_TYPE x_shift = vec - maxVal * VEC_TYPE::Ones(dim);

	double lse_x = log(x_shift.array().exp().sum()) + maxVal;
	
	return lse_x;
  } // end of lse() function


VEC_TYPE VarParam::lse_rows(MAT_TYPE X, int N_rows, int D_cols) {

	VEC_TYPE lse_vec = VEC_TYPE::Zero(N_rows);
	for (int n = 0; n < N_rows; n++) {
		lse_vec(n) = lse(X.row(n), D_cols);
	}
	return lse_vec;

} // end of lse_rows() function



void VarParam::nextIter() {
	this->curr++;
} // end of nextIter() function


bool VarParam::elboConverge() {

	int i = this->curr;

	double diff = this->L(i) - this->L(i - 1);

	if (i % 20 == 0) {
		cout << "It:\t" << i;
		std::cout << std::setprecision(4) << std::fixed;
		cout << "\tLB:\t" << this->L(i) << "\tdelta:\t" << diff << endl;
	}
	
	bool converge = abs(diff) < this->tol;

	return converge;

} // end elboConverge() function



void VarParam::eStep() {

	int n, k;
	int N = getN();
	int K = getClusters();

	VEC_TYPE xVx = VEC_TYPE::Zero(K);

	VEC_TYPE psi_a = digamma(this->a_k.array());
	VEC_TYPE psi_b = digamma(this->b_k.array());

	VEC_TYPE ONES_K = VEC_TYPE::Ones(K);         // K-dim vector of ones
	// VEC_TYPE ONES_N = VEC_TYPE::Ones(N);         // N-dim vector of ones


	MAT_TYPE X_mu   = (this->X) * (this->mu_k);  // (N x K) : X * mu_k

	for (n = 0; n < N; n++) {
		
		xVx.fill(0);
		VEC_TYPE x_n   = this->X.row(n);
		VEC_TYPE xmu_n = X_mu.row(n);

		for (k = 0; k < K; k++) {

			xVx(k) = (x_n.transpose() * (*this->Vk_inv_it) * x_n).value();

			advance(this->Vk_inv_it, 1);

		} // end inner for
		this->Vk_inv_it = this->V_k_inv.begin(); // bring back to front of list

		// compute n-th row of log_rho_nk (K x 1)
		VEC_TYPE t1 = xmu_n - ((this->alpha(n) + this->phi(n)) * ONES_K); 
		VEC_TYPE t2 = log(2 * M_PI) * ONES_K - psi_a + psi_b;
		VEC_TYPE tau_prod = this->tau_k.array().cwiseProduct(((this->y(n) * 
			ONES_K).array() - (x_n.transpose() * this->m_k).transpose().array()).square());

		this->log_r_nk.row(n) = t1 - 0.5 * (t2 + xVx + tau_prod);

	} // end outer for

	/*  compute r_nk = divide each element of rho_nk by the sum of the
		corresponding row; each row consists of 'responsibilities'       */
    
	VEC_TYPE logZ  = lse_rows(this->log_r_nk, N, K);
	this->log_r_nk = (this->log_r_nk).colwise() - logZ; // normalize log_r_nk

	this->r_nk     = (this->log_r_nk).array().exp();    // exp() to recover r_nk
	
	this->N_k = this->r_nk.colwise().sum();



} // end eStep() function



void VarParam::mStep() {

	// printf("running m-step.\n");

	int n, k;

	int N = getN(); 
	int D = getCovDim();
	int K = getClusters();

	VEC_TYPE xQx(K);                             // K-dim vector used in xi
	VEC_TYPE ONES_K = VEC_TYPE::Ones(K);         // K-dim vector of ones
	VEC_TYPE ONES_N = VEC_TYPE::Ones(N);         // N-dim vector of ones
	MAT_TYPE I_D    = MAT_TYPE::Identity(D, D);  // (D x D) identity matrix
 
	MAT_TYPE X_mu   = (this->X) * (this->mu_k);  // (N x K) : X * mu_k

	MAT_TYPE rl_nk_xx = MAT_TYPE::Zero(D, D);    // used in q(gamma)
	MAT_TYPE r_nk_xx  = MAT_TYPE::Zero(D, D); 	 // used in q(beta|tau)
	MAT_TYPE r_x      = MAT_TYPE::Zero(D, D);    // used in rl_nk_xx, r_nk_xx

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
	  *     (1.3) q(tau)       DONE
	  * 
	  * (2) update posterior means
	  *     (2.1) beta_k       DONE
	  *     (2.2) tau_k        DONE
	  *     (2.3) gamma_k      DONE
	  ----------------------------------------------------------------------  */

	/* q(gamma) : Q_k, Q_k^{-1}, eta_k, mu_k -------------------------------- */
	for (k = 0; k < K; k++) {

		rl_nk_xx = MAT_TYPE::Zero(D, D); // reset entries to 0 every iter
		r_nk_xx  = MAT_TYPE::Zero(D, D); // reset entries to 0 every iter

		// calculate intermediate quantities for q(gamma), q(beta|tau)
		for (n = 0; n < N; n++) {
			// VEC_TYPE x_n = this->X.row(n);
			// x_n * x_n' can be precomputed 
			r_x      = this->r_nk(n,k) * ((*this->xn_mat_it)).array();
			rl_nk_xx = rl_nk_xx.array() + this->lambda(n, k) * r_x.array(); 
			r_nk_xx  = r_nk_xx + r_x;

			advance(this->xn_mat_it, 1);
		} // end inner for

		this->xn_mat_it = this->xn_mat.begin(); 

		/* (1.1) q(gamma)    : Q_k, Q_k_inv, eta_k, mu_k -------------------- */
		
		VEC_TYPE rnk_k = this->r_nk.col(k);                         // (N x 1)
		VEC_TYPE lam_k = this->lambda.col(k);                       // (N x 1)
 		VEC_TYPE alpha_lam_k = 0.5 * ONES_N + 
 								2 * lam_k.cwiseProduct(this->alpha);  

		*this->Qk_it      = I_D.array() + 2 * rl_nk_xx.array();
		*this->Qk_inv_it  = (*this->Qk_it).ldlt().solve(I_D);
		
		this->eta_k.col(k) = this->X.transpose() * 
								(rnk_k.cwiseProduct(alpha_lam_k));   // (D x 1)
		this->mu_k.col(k)  = (*this->Qk_inv_it) * this->eta_k.col(k);

		/* (1.2) q(beta|tau) : V_k, V_k_inv, zeta_k, m_k -------------------- */
		VEC_TYPE rnk_y = rnk_k.cwiseProduct(this->y);

		*this->Vk_it     = this->Lambda_0 + r_nk_xx;
		*this->Vk_inv_it = (*this->Vk_it).ldlt().solve(I_D);
		
		this->zeta_k.col(k) = this->Lambda0_m0.array() + 
									(this->X.transpose() * rnk_y).array();
		this->m_k.col(k)    = (*this->Vk_inv_it) * this->zeta_k.col(k);


		// update b_k
		this->b_k(k) = - ((this->zeta_k.col(k)).transpose() * 
			((*this->Vk_inv_it) * (this->zeta_k.col(k)))).value() + 
					(rnk_k.cwiseProduct(this->y2)).sum();


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


	/* (1.3) update q(tau)   : a_k, b_k ------------------------------------- */
	
	this->a_k = (this->a_0 * ONES_K).array() + (0.5 * this->N_k).array();
	/*
	for (k = 0; k < K; k++) {
		VEC_TYPE rnk_k = this->r_nk.col(k);                         // (N x 1)
		this->b_k(k) = - ((this->zeta_k.col(k)).transpose() * 
			((*this->Vk_inv_it) * (this->zeta_k.col(k)))).value() + 
					(rnk_k.cwiseProduct(this->y2)).sum();
		advance(this->Vk_inv_it, 1);
	}
	*/

	// this->Vk_inv_it = V_k_inv.begin(); 

	this->b_k = this->b_0.array() + 
				(0.5 * (this->b_k + this->m0_Lambda0_m0 * ONES_K)).array();

	
	/** IN PROGRESS -------------------------------------------------------- **/

	/* (2.1) update beta_k  ------------------------------------------------- */
	this->beta_k = this->m_k;

	/* (2.2) update tau_k   ------------------------------------------------- */
	this->tau_k = (this->a_k).cwiseQuotient(this->b_k);

	/* (2.3) update gamma_k ------------------------------------------------- */
	this->gamma_k = this->mu_k;

	// update current iteration
	//this->curr++;


} // end mStep() function


void VarParam::elbo() {

	/*  Calculate the 7 following expectations whose sum is the ELBO for the
		the current iteration fo the CAVI algorithm

	    (1) E [ ln p(y | X, beta, tau, Z) ] 
	    (2) E [ ln p(Z | X, gamma) ] 
	    (3) E [ ln p(gamma) ] 
		(4) E [ ln p(beta, tau) ] 
	    
		(5) E [ ln q(Z) ]
		(6) E [ ln q(beta, tau) ]
		(7) E [ ln q(gamma) ]
	
	*/



	int n, k;
	int N = getN();
	int D = getCovDim();
	int K = getClusters();

	VEC_TYPE xVx    = VEC_TYPE::Zero(K);     // (K x 1)
	VEC_TYPE ONES_K = VEC_TYPE::Ones(K);     // K-dim vector of ones

	// initialize the 7 expectations that are computed in the ELBO
	double e_ln_p_y = 0, e_ln_p_z = 0, e_ln_p_gamma = 0, e_ln_p_beta_tau = 0;
	double e_ln_q_z = 0, e_ln_q_beta_tau = 0, e_ln_q_gamma = 0;

	// commonly computed quantities
	VEC_TYPE psi_a = digamma(this->a_k.array());    // (K x 1)
	VEC_TYPE psi_b = digamma(this->b_k.array());    // (K x 1)
	
	MAT_TYPE X_mu  = (this->X) * (this->mu_k);      // (N x K) : X * mu_k

	MAT_TYPE diff  = (this->m_k).colwise() - this->m_0; // (N x K), for (4)
	
	
	VEC_TYPE e1 = VEC_TYPE::Zero(N);
	VEC_TYPE e2 = VEC_TYPE::Zero(K);
	
	VEC_TYPE e4_diff  = VEC_TYPE(K); // (4) store quadratic term
	VEC_TYPE e4_trace = VEC_TYPE(K); // (4) store trace term

	VEC_TYPE e6 = VEC_TYPE(K); // (6) store log(det(V_k_inv)) terms

	VEC_TYPE e7 = VEC_TYPE(K); // (7) store log(det(Q_k))


	// calculate (1)
	for (n = 0; n < N; n++) {

		xVx.fill(0);
		VEC_TYPE x_n   = this->X.row(n);
		VEC_TYPE xmu_n = X_mu.row(n);
		VEC_TYPE rnk_n = this->r_nk.row(n);

		for (k = 0; k < K; k++) {
			
			// calculation for (1) 
			xVx(k) = (x_n.transpose() * (*this->Vk_inv_it) * x_n).value();
			advance(this->Vk_inv_it, 1);

		} // end inner for (k)

		this->Vk_inv_it = this->V_k_inv.begin(); // reset V_k_inv iterator

		VEC_TYPE t2 = log(2 * M_PI) * ONES_K - psi_a + psi_b + xVx;
		VEC_TYPE tau_prod = this->tau_k.array().cwiseProduct(((this->y(n) * 
			ONES_K).array() - (x_n.transpose() * this->m_k).transpose().array()).square());

		e1(n) = (rnk_n.array().cwiseProduct((t2 + tau_prod).array())).sum();

	} // end outer for (n)


	// calculate (2), (4), (6), (7)
	for (k = 0; k < K; k++) {

		// calculation for (2)
		VEC_TYPE rnk_k = this->r_nk.col(k);  // (N x 1)
		VEC_TYPE xmu_k = X_mu.col(k);        // (N x 1)
		e2(k) = (rnk_k.array().cwiseProduct((xmu_k - this->alpha - 
												this->phi).array())).sum();


		// calculation for (4) -- update e4_diff, e4_trace
		VEC_TYPE diff_k = diff.col(k);
		e4_diff(k)  = (diff_k.transpose() * 
								(this->Lambda_0) * diff_k).value();
		e4_trace(k) = (this->Lambda_0 * (*this->Vk_inv_it)).trace();

		// calculation for (6) -- update e6
		e6(k) = log((*this->Vk_it).determinant());


		// calculation for (7) -- update e7
		e7(k) = log((*this->Qk_it).determinant());

		advance(this->Vk_inv_it, 1);
		advance(this->Qk_it, 1);
		advance(this->Vk_it, 1);
		

	} // end for (k)

	this->Qk_it     = this->Q_k.begin();     // reset Q_k iterator
	this->Vk_it     = this->V_k.begin();     // reset V_k iterator
	this->Vk_inv_it = this->V_k_inv.begin(); // reset V_k_inv iterator



	// calculation for (4)
	double e3_1 = (this->a_0 + (0.5 * D - 1) * 
					ONES_K).array().cwiseProduct((psi_a - psi_b).array()).sum();

	double e3_2 = - K * (0.5 * D * log(2 * M_PI) - this->lgd_Lambda_0 - 
				  		this->a_0(0) * this->log_b0(0) + this->lg_a0(0));

	e4_diff = (this->tau_k).array().cwiseProduct((e4_diff + this->b_0).array());

	VEC_TYPE e6_const = ((0.5 * D - 1) * ONES_K + 
			this->a_k).array().cwiseProduct((psi_a - psi_b).array()).array() + 
			(this->a_k).array().cwiseProduct(this->b_k.array().log()) - 
			(this->a_k).array() - (this->a_k).array().lgamma();

	e_ln_p_y        = -0.5 * e1.array().sum();
	e_ln_p_z        = e2.array().sum();
	e_ln_p_gamma    = -0.5 * K * D * log(2 * M_PI) - 
						  0.5 * ((this->mu_k).array().square().sum());
	e_ln_p_beta_tau = e3_1 + e3_2 - 0.5 * (e4_diff + e4_trace).array().sum();
	e_ln_q_z  	    = (this->r_nk.array().cwiseProduct(
										this->log_r_nk.array())).array().sum();
	e_ln_q_beta_tau = (e6_const + e6).array().sum() - 
						0.5 * K * D * (log(2 * M_PI) + 1);

	e_ln_q_gamma    = -0.5 * K * D * (log(2 * M_PI) + 1) + e7.array().sum();					

	/*
	std::cout.precision(8);
	std::cout << "e_ln_p_y        = " << e_ln_p_y << endl;
	std::cout << "e_ln_p_z        = " << e_ln_p_z << endl;
	std::cout << "e_ln_p_gamma    = " << e_ln_p_gamma << endl;
	std::cout << "e_ln_p_beta_tau = " << e_ln_p_beta_tau << endl;
	std::cout << "e_ln_q_z        = " << e_ln_q_z << endl;
	std::cout << "e_ln_q_beta_tau = " << e_ln_q_beta_tau << endl;
	std::cout << "e_ln_q_gamma    = " << e_ln_q_gamma << endl;
	*/

	// update the ELBO
	this->L(this->curr) = e_ln_p_y + e_ln_p_z + e_ln_p_gamma + e_ln_p_beta_tau - 
									e_ln_q_z - e_ln_q_gamma - e_ln_q_beta_tau;

} // end of elbo() function



// end of VarParam.cpp file

int main(void) {
	printf("inside variational parameter class file\n");
	return 0;
}

