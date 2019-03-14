


// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>
#include <numeric>

using namespace std;
using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;

#define VEC_TYPE  VectorXd
#define MAT_TYPE  MatrixXd

#define MAP_MAT Eigen::Map<Eigen::MatrixXd> 
#define MAP_VEC Eigen::Map<Eigen::VectorXd> 


// [[Rcpp::export]]
SEXP mat_list_ops(int K, int D, VEC_TYPE x) {     		 // (1 x 1)

	int k;
	MAT_TYPE I_D = MAT_TYPE::Identity(D, D);

	// MAT_TYPE I_D = MAT_TYPE::Zero(D, D);
	// I_D.setZero(D, D);

	
	list<MAT_TYPE> Q_k;
	list<MAT_TYPE> prod;

	// populate the list with matrices
	for (k = 0; k < K; k++) {
		MAT_TYPE mat_k = (k + 1) * I_D;
		Q_k.push_back(mat_k);
	}

	list<MAT_TYPE>::iterator it = Q_k.begin();
	// advance(it, 2); // Advance the iterator by 2 positions,

	// iterate the list of matrices and multiply each matrix with vector x
	for (k = 0; k < K; k++) {
		VEC_TYPE y = (*it) * x;
		prod.push_back(y);
		advance(it, 1);
	}


	return wrap(zero_vec);
	
} // end mat_list_ops() function

