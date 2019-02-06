
## betaSimAnalysis.R


## read in results of betaSimulations.R:
##     (1) overlayed densities (true, vb1, vb2)
##     TODO: (2) point estimates (true - approx)
##     (3) convergence related figures
##         (i)   elbo plots
##         (ii)  iterative changes in approximating density

## first set directory to /code directory

source("globals.R")

setwd(HOME_DIR)

source(DENSITY)
source(BETA_1)
source(BETA_2)
source(BETA_2_MODIFIED)
source(paste(COV_INDEP, MISC_FUNCS, sep = '/')) 

# read in data files for each of the simulations -------------------------------

# (1) X : (200 x 2)
theta1_1 = readRDS(file = "beta_results/theta_indep_beta1_200.rds")
theta1_2 = readRDS(file = "beta_results/theta_dep_beta1_200.rds")


# (2) X : (500 x 4)
theta2_1 = readRDS(file = "beta_results/theta_indep_beta2_500.rds")
theta2_2 = readRDS(file = "beta_results/theta_dep_beta2_500.rds")


# (3) X : (1000 x 10)
theta3_1 = readRDS(file = "beta_results/theta_indep_beta2_mod_1000.rds")
theta3_2 = readRDS(file = "beta_results/theta_dep_beta2_mod_1000.rds")


# ------------------------------------------------------------------------------

y_grid    = seq(0, 1, len = 500)            # grid to feed into densities
approx_d  = list(py_0, py_bouch)            # list of approx density functions
den_label = c("cov-indep", "cov-dep (b)")   # labels for each approx density
K         = 4                               # number of clusters


# overlay densities for simulation ---------------------------------------------


# (1) --------------------------------------------------------------------------
N           = 200                           # number of observations
synth_data1 = old_beta1(N, seed = 123)
y1          = synth_data1$y                 # (N x 1)
X1          = synth_data1$X                 # (N x 2)
param_mat1  = synth_data1$shape             # (N x 2)
theta1      = list(theta1_1, theta1_2)      # list of var. params for each alg

n_set1   = c(10, 30, 42, 88)
p_list1 = queryDensity(X1, n_set1, param_mat1, dbeta,
                      approx_d, den_label, theta1, K)

multiplot(plotlist = p_list1, cols = 2)

# (2) --------------------------------------------------------------------------
N           = 500                                # number of observations
K           = 4                                  # number of clusters
synth_data2 = beta2(N)
y2          = synth_data2$y                      # (N x 1)
X2          = synth_data2$X                      # (N x 2)
param_mat2  = synth_data2$shape                  # (N x 2)
theta2      = list(theta2_1, theta2_2) 

n_set2  = c(100, 200, 312, 400)
p_list2 = queryDensity(X2, n_set2, param_mat2, dbeta,
                      approx_d, den_label, theta2, K)

multiplot(plotlist = p_list2, cols = 2)


# (3) --------------------------------------------------------------------------
N           = 1000                               # number of observations
K           = 4                                  # number of clusters
synth_data3 = beta2_modified(N)
y3          = synth_data3$y                      # (N x 1)
X3          = synth_data3$X                      # (N x 2)
param_mat3  = synth_data3$shape                  # (N x 2)
theta3      = list(theta3_1, theta3_2)

n_set3  = c(250, 500, 111, 291)
p_list3 = queryDensity(X3, n_set3, param_mat3, dbeta,
                      approx_d, den_label, theta3, K)

multiplot(plotlist = p_list3, cols = 2)

plotELBO(theta3_1)
plotELBO(theta3_2)






# end of bestaSimAnalysis.R
