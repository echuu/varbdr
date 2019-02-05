
## simulations.R

HOME_DIR = "~/varbdr/code"              # linux
HOME_DIR = "C:/Users/chuu/varbdr/code"  # windows

setwd(HOME_DIR)
source("beta1.R")    # beta1() used to generate synthetic data

## run the VB algorithms and compare performance via:
##     (1) overlayed densities (true, vb1, vb2)
##     (2) point estimates (true - approx)
##     (3) convergence related figures
##         (i)   elbo plots
##         (ii)  iterative changes in approximating density


# (0) simulate first set of synthetic data (beta1)
N           = 1000                     # number of observations
K           = 5                        # number of clusters
synth_data1 = beta1(N)
y1          = synth_data1$y            # (N x 1)
X1          = synth_data1$X            # (N x 2)
param_mat   = synth_data1$shape        # (N x 2)

# (1) covariate-INDEPENENT vb
setwd(HOME_DIR)                        # navigate back to parent directory
source("cov_indep/varbdr.R")           # navigate to directory that has algo (1)
theta1 = varbdr(y = y1, X = X1, K = K) # store results of CAVI converges

# (2) covariate-DEPENDENT vb
setwd(HOME_DIR)                         # navigate back to parent directory
source("cov_dep/varbdr.R")             # navigate to directory that has algo (1)
theta2 = varbdr(y = y1, X = X1, K = K) # store results of CAVI converges


# (3) overlay corresponding densities obtained from steps (1) and (2)
setwd(HOME_DIR)                 # navigate back to parent directory
source("density.R")

n = 421
y_grid = seq(0, 1, len = 500)
params = list(shape1 = param_mat[n,1], shape2 = param_mat[n,2])

p1 = plotDensities(y_grid, X1[n,], dbeta, params, py_0, theta1, K)
p1$overlay

p2 = plotDensities(y_grid, X1[n,], dbeta, params, py_bouch, theta2, K)
p2$overlay


# ------------------------------------------------------------------------------
library(reshape2)
source("density.R")

n = 421
x       = X1[n,]
y_grid  = seq(0, 1, len = 500)
params = list(shape1 = param_mat[n,1], shape2 = param_mat[n,2])
theta     = list(theta1, theta2)          # list of var. params for each alg
approx_d  = list(py_0, py_bouch)          # list of approx density functions
den_label = c("cov-indep", "cov-dep (b)") # labels for each approx density

compareDensities(y_grid, x, 
                 dbeta, params,
                 approx_d, den_label, 
                 theta, K)




# (4) look at convergence of elbo for (1) and (2)

setwd(HOME_DIR)                        # navigate back to parent directory
source("cov_indep/misc.R") 

plotELBO(theta1)

plotELBO(theta2)



# end of simulations.R
