
## betaSimulations.R

## first set directory to /code directory

source("globals.R")

setwd(HOME_DIR)
source(BETA_1)                                   # (1) beta1 data sim
source(BETA_2)                                   # (2) beta2 data sim
source(BETA_2_MODIFIED)                          # (3) beta2 (mod) data sim

## run the VB algorithms for each of the data generating process


# (1) beta1 simulation ---------------------------------------------------------
N           = 200                                # number of observations
K           = 4                                  # number of clusters
synth_data1 = old_beta1(N, seed = 123)
y1          = synth_data1$y                      # (N x 1)
X1          = synth_data1$X                      # (N x 2)
param_mat1  = synth_data1$shape                  # (N x 2)

# (1.1) covariate-INDEPENENT vb
source(paste(COV_INDEP, VARBDR, sep = '/'))      # load cov-indep varbdr.R
theta1_1 = varbdr(y = y1, X = X1, K = K)           # store CAVI results

# (1.2) covariate-DEPENDENT vb
source(paste(COV_DEP, VARBDR, sep = '/'))        # load cov-dep varbdr.R
theta1_2 = varbdr(y = y1, X = X1, K = K)         # store CAVI results


# (2) beta2 simulation ---------------------------------------------------------
N           = 500                                # number of observations
K           = 4                                  # number of clusters
synth_data2 = beta2(N)
y2          = synth_data2$y                      # (N x 1)
X2          = synth_data2$X                      # (N x 2)
param_mat2  = synth_data2$shape                  # (N x 2)

# (2.1) covariate-INDEPENENT vb
source(paste(COV_INDEP, VARBDR, sep = '/'))      # load cov-indep varbdr.R
theta2_1 = varbdr(y = y2, X = X2, K = K)         # store CAVI results

# (2.2) covariate-DEPENDENT vb
source(paste(COV_DEP, VARBDR, sep = '/'))        # load cov-dep varbdr.R
theta2_2 = varbdr(y = y2, X = X2, K = K)         # store CAVI results


# (3) beta2 (MODIFIED) simulation ----------------------------------------------
N           = 1000                               # number of observations
K           = 4                                  # number of clusters
synth_data3 = beta2_modified(N)
y3          = synth_data3$y                      # (N x 1)
X3          = synth_data3$X                      # (N x 2)
param_mat3  = synth_data3$shape                  # (N x 2)

# (3.1) covariate-INDEPENENT vb
source(paste(COV_INDEP, VARBDR, sep = '/'))      # load cov-indep varbdr.R
theta3_1 = varbdr(y = y3, X = X3, K = K)         # store CAVI results

# (3.2) covariate-DEPENDENT vb
source(paste(COV_DEP, VARBDR, sep = '/'))        # load cov-dep varbdr.R
theta3_2 = varbdr(y = y3, X = X3, K = K)         # store CAVI results



## save variational parameters as R objects ------------------------------------

saveRDS(theta1_1, file = "beta_results/theta_indep_beta1_200.RDS")
saveRDS(theta1_2, file = "beta_results/theta_dep_beta1_200.RDS")

saveRDS(theta2_1, file = "beta_results/theta_indep_beta2_500.RDS")
saveRDS(theta2_2, file = "beta_results/theta_dep_beta2_500.RDS")

saveRDS(theta3_1, file = "beta_results/theta_indep_beta2_mod_1000.RDS")
saveRDS(theta3_2, file = "beta_results/theta_dep_beta2_mod_1000.RDS")


# end of betaSimulations.R
