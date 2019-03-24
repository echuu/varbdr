


// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]

// #include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>
#include <numeric>
//#include <Rmath.h>

#include <unsupported/Eigen/SpecialFunctions>

#include<math.h>


//#include <boost/math/special_functions/digamma.hpp>
// #include "VarParam.h"

using namespace std;
using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;

#define VEC_TYPE  VectorXd
#define MAT_TYPE  MatrixXd

#define MAP_MAT Eigen::Map<Eigen::MatrixXd> 
#define MAP_VEC Eigen::Map<Eigen::VectorXd> 




// [[Rcpp::export]]
double lse(VEC_TYPE vec, int dim) {

	


// [[Rcpp::export]]
SEXP lse_rows(MAT_TYPE X, int N, int D) { 		 // (1 x 1)
 
 	VEC_TYPE lse_vec = VEC_TYPE::Zero(N);
	for (int n = 0; n < N; n++) {
		lse_vec(n) = lse(X.row(n), D);
	}
	
	return(wrap(lse_vec));
} // end mat_list_ops() function

