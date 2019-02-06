
## simulations.R

source("globals.R")

setwd(HOME_DIR)
source(BETA_1)                                   # data generating scheme # 1

## run the VB algorithms and compare performance via:
##     (1) overlayed densities (true, vb1, vb2)
##     (2) point estimates (true - approx)
##     (3) convergence related figures
##         (i)   elbo plots
##         (ii)  iterative changes in approximating density


# (0) simulate first set of synthetic data (beta1)
N           = 200                                # number of observations
K           = 5                                  # number of clusters
synth_data1 = old_beta1(N, seed = 123)
y1          = synth_data1$y                      # (N x 1)
X1          = synth_data1$X                      # (N x 2)
param_mat1  = synth_data1$shape                  # (N x 2)

# (1) covariate-INDEPENENT vb
source(paste(COV_INDEP, VARBDR, sep = '/'))      # load cov-indep varbdr.R
theta1_1 = varbdr(y = y1, X = X1, K = K)           # store CAVI results

# (2) covariate-DEPENDENT vb
source(paste(COV_DEP, VARBDR, sep = '/'))        # load cov-dep varbdr.R
theta1_2 = varbdr(y = y1, X = X1, K = K)         # store CAVI results


# (3) overlay corresponding densities obtained from steps (1) and (2)
source(DENSITY)

n = 30
y_grid = seq(0, 1, len = 500)
params1 = list(shape1 = param_mat1[n,1], shape2 = param_mat1[n,2])

p1 = plotDensities(y_grid, X1[n,], dbeta, params1, py_0, theta1_1, K)
p1$overlay

p2 = plotDensities(y_grid, X1[n,], dbeta, params1, py_bouch, theta1_2, K)
p2$overlay


# ------------------------------------------------------------------------------

# same for each beta algorithm
y_grid  = seq(0, 1, len = 500)
approx_d  = list(py_0, py_bouch)          # list of approx density functions
den_label = c("cov-indep", "cov-dep (b)") # labels for each approx density


n = 10
x       = X1[n,]
params1 = list(shape1 = param_mat1[n,1], shape2 = param_mat1[n,2])
theta1  = list(theta1_1, theta1_2)          # list of var. params for each alg

p_10 = compareDensities(y_grid, x, 
                        dbeta, params1,
                        approx_d, den_label, 
                        theta1, K)


n = 30
x       = X1[n,]
params1 = list(shape1 = param_mat1[n,1], shape2 = param_mat1[n,2])
theta1  = list(theta1_1, theta1_2)          # list of var. params for each alg

p_30 = compareDensities(y_grid, x, 
                        dbeta, params1,
                        approx_d, den_label, 
                        theta1, K)


n = 42
x       = X1[n,]
params1 = list(shape1 = param_mat1[n,1], shape2 = param_mat1[n,2])
theta1  = list(theta1_1, theta1_2)          # list of var. params for each alg

(p_42 = compareDensities(y_grid, x, 
                        dbeta, params1,
                        approx_d, den_label, 
                        theta1, K))


n = 88
x       = X1[n,]
params1 = list(shape1 = param_mat1[n,1], shape2 = param_mat1[n,2])
theta1  = list(theta1_1, theta1_2)          # list of var. params for each alg

(p_88 = compareDensities(y_grid, x, 
                         dbeta, params1,
                         approx_d, den_label, 
                         theta1, K))


multiplot(p_10, p_30, p_42, p_88, cols = 2)


# (4) look at convergence of elbo for (1) and (2)

setwd(HOME_DIR)                        # navigate back to parent directory
source("cov_indep/misc.R") 

plotELBO(theta1_1)
plotELBO(theta1_2)

## beta2 simulations -----------------------------------------------------------
source(BETA_2)
N           = 500                                # number of observations
K           = 4                                  # number of clusters
synth_data2 = beta2(N)
y2          = synth_data2$y                      # (N x 1)
X2          = synth_data2$X                      # (N x 2)
param_mat2  = synth_data2$shape                  # (N x 2)

# (1) covariate-INDEPENENT vb
source(paste(COV_INDEP, VARBDR, sep = '/'))      # load cov-indep varbdr.R
theta2_1 = varbdr(y = y2, X = X2, K = K)         # store CAVI results

# (2) covariate-DEPENDENT vb
source(paste(COV_DEP, VARBDR, sep = '/'))        # load cov-dep varbdr.R
theta2_2 = varbdr(y = y2, X = X2, K = K)         # store CAVI results

n = 100
x       = X2[n,]
params2 = list(shape1 = param_mat2[n,1], shape2 = param_mat2[n,2])
theta2  = list(theta2_1, theta2_2)               # list of var. params

compareDensities(y_grid, x, 
                 dbeta, params2,
                 approx_d, den_label, 
                 theta2, K)

plotELBO(theta2_1)
plotELBO(theta2_2)

# save R objects for N = 500, K = 4
saveRDS(theta2_1, file = "theta_indep_beta2_500.RDS")
saveRDS(theta2_2, file = "theta_dep_beta2_500.RDS")


theta2_1 = readRDS(file = "theta_indep_beta2_500.rds")
theta2_2 = readRDS(file = "theta_dep_beta2_500.rds")


n = 100
x       = X2[n,]
params2 = list(shape1 = param_mat2[n,1], shape2 = param_mat2[n,2])
theta2  = list(theta2_1, theta2_2)               # list of var. params

p_100 = compareDensities(y_grid, x, 
                         dbeta, params2,
                         approx_d, den_label, 
                         theta2, K)

n = 200
x       = X2[n,]
params2 = list(shape1 = param_mat2[n,1], shape2 = param_mat2[n,2])
theta2  = list(theta2_1, theta2_2)               # list of var. params

p_200 = compareDensities(y_grid, x, 
                         dbeta, params2,
                         approx_d, den_label, 
                         theta2, K)


n = 312
x       = X2[n,]
params2 = list(shape1 = param_mat2[n,1], shape2 = param_mat2[n,2])
theta2  = list(theta2_1, theta2_2)               # list of var. params

(p_312 = compareDensities(y_grid, x, 
                         dbeta, params2,
                         approx_d, den_label, 
                         theta2, K))

n = 400
x       = X2[n,]
params2 = list(shape1 = param_mat2[n,1], shape2 = param_mat2[n,2])
theta2  = list(theta2_1, theta2_2)               # list of var. params

p_400 = compareDensities(y_grid, x, 
                         dbeta, params2,
                         approx_d, den_label, 
                         theta2, K)


multiplot(p_100, p_200, p_312, p_400, cols = 2)


# ------------------------------------------------------------------------------
 

source(BETA_2_MODIFIED)
N           = 1000                               # number of observations
K           = 4                                  # number of clusters
synth_data3 = beta2_modified(N)
y3          = synth_data3$y                      # (N x 1)
X3          = synth_data3$X                      # (N x 2)
param_mat3  = synth_data3$shape                  # (N x 2)

# (1) covariate-INDEPENENT vb
source(paste(COV_INDEP, VARBDR, sep = '/'))      # load cov-indep varbdr.R
theta3_1 = varbdr(y = y3, X = X3, K = K)         # store CAVI results

# (2) covariate-DEPENDENT vb
source(paste(COV_DEP, VARBDR, sep = '/'))        # load cov-dep varbdr.R
theta3_2 = varbdr(y = y3, X = X3, K = K)         # store CAVI results


saveRDS(theta3_1, file = "theta_indep_beta2_mod_1000.RDS")
saveRDS(theta3_2, file = "theta_dep_beta2_mod_1000.RDS")


n = 250
x       = X3[n,]
params3 = list(shape1 = param_mat3[n,1], shape2 = param_mat3[n,2])
theta3  = list(theta3_1, theta3_2)               # list of var. params

(p_250 = compareDensities(y_grid, x, 
                         dbeta, params3,
                         approx_d, den_label, 
                         theta3, K))

n = 500
x       = X3[n,]
params3 = list(shape1 = param_mat3[n,1], shape2 = param_mat3[n,2])
theta3  = list(theta3_1, theta3_2)               # list of var. params

(p_500 = compareDensities(y_grid, x, 
                         dbeta, params3,
                         approx_d, den_label, 
                         theta3, K))

n = 111
x       = X3[n,]
params3 = list(shape1 = param_mat3[n,1], shape2 = param_mat3[n,2])
theta3  = list(theta3_1, theta3_2)               # list of var. params

(p_111 = compareDensities(y_grid, x, 
                         dbeta, params3,
                         approx_d, den_label, 
                         theta3, K))

n = 291
x       = X3[n,]
params3 = list(shape1 = param_mat3[n,1], shape2 = param_mat3[n,2])
theta3  = list(theta3_1, theta3_2)               # list of var. params

(p_291 = compareDensities(y_grid, x, 
                         dbeta, params3,
                         approx_d, den_label, 
                         theta3, K))

multiplot(p_250, p_500, p_111, p_291, cols = 2)



plotELBO(theta3_1)
plotELBO(theta3_2)



# end of simulations.R
