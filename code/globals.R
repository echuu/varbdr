
## globals.R

# libraries used
library(ggplot2)
library(reshape2)
library(matrixcalc)
library(dplyr)

# libraries needed to run the C++ code
#library("Rcpp")
library("microbenchmark")
#library("RcppEigen")
#library("RcppArmadillo")
#library("RcppParallel")


# home directory of all the code
HOME_DIR = "~/varbdr/code"                # linux
# HOME_DIR = "C:/Users/chuu/varbdr/code"      # windows

# directories of each of the different algorithms
COV_INDEP = paste(HOME_DIR, "/cov_indep", sep = '')
COV_DEP   = paste(HOME_DIR, "/cov_dep", sep = '')


# C++ code directories
CPP_DIR   = paste(HOME_DIR, "cpp_code", sep = '')


# cpp code for faster matrix operations
VB_CPP = "vb_calcs.cpp"


# file names (mostly common to all algorithms)
INIT_PRIORS      =  "initPriors.R"
INIT_VAR_PARAMS  =  "initVarParams.R"
VARBDR           =  "varbdr.R"
FAST_VARBDR      =  "fast_varbdr.R"
E_STEP           =  "eStep.R"
M_STEP           =  "mStep.R"
FAST_M_STEP      =  "fast_mStep.R"
ELBO             =  "elbo.R"
MISC_FUNCS       =  "misc.R"
DENSITY          =  "density.R"


# ID's for data generating schemes
BETA     = 0
DP_MIX1  = 1
DP_MIX2  = 2
BIMODAL1 = 3

# data-generating schemes
BETA_1           = "data_gen/beta1.R"
BETA_2           = "data_gen/beta2.R"
BETA_2_MODIFIED  = "data_gen/beta2_modified.R"
DP_BDR           = "data_gen/dp_bdr.R"
UNI_CDE          = 'data_gen/uni_cde.R'    



