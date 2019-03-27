


// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]

// #include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>
#include <numeric>
#include "ParallelOps.h"

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



/// Basically replacing:
// [[Rcpp::export]]
void sequential_for (MAT_TYPE X) {

	MAT_TYPE I_D = MAT_TYPE::Identity(X.rows(), X.rows());

    for (int i = 0; i < 10; ++i)
        MAT_TYPE X_inv = X.llt().solve(I_D);
}

/// By:
// [[Rcpp::export]]
void threaded_for (MAT_TYPE X) {

	//int start = 0;
	//int end = 10;
	MAT_TYPE I_D = MAT_TYPE::Identity(X.rows(), X.rows());

    parallel_for(10, [&](int start, int end){ 
        for(int i = start; i < end; ++i)
            MAT_TYPE X_inv = X.llt().solve(I_D);
    } );
}


// [[Rcpp::export]]
MAT_TYPE test (MAT_TYPE X, MAT_TYPE mu_k) {

		VEC_TYPE x_n = X.row(0);

		MAT_TYPE M = (x_n.transpose() * mu_k).transpose();

		return M;

}

// [[Rcpp::export]]
SEXP fastInverse(MAT_TYPE X) {

	MAT_TYPE I_D = MAT_TYPE::Identity(X.rows(), X.rows());

	MAT_TYPE X_inv = X.llt().solve(I_D);

	return wrap(X_inv);

}

// [[Rcpp::export]]
SEXP slowInverse(MAT_TYPE X) {

	// MAT_TYPE I_D = MAT_TYPE::Identity(X.rows(), X.rows());

	MAT_TYPE X_inv = X.inverse();

	return wrap(X_inv);

}


