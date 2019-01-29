
## testFunctions.R -- test/troubleshoot functions (covariate-DEPENDENT)


# generate data ----------------------------------------------------------------

N        = 100  # number of observations
D        = 2    # number of covariates 
K        = 5    # number of clusters
max_iter = 200  # max number of iterations 

# generate X1,...,Xp ~ Unif[0.05, 0.95]
## initialize storage for response, covariates
y = numeric(N)       # (N x 1) vector of response variables
X = matrix(0, N, D)  # (N x D) matrix; covariates for y_n are stored row-wise  

shape_mat = matrix(0, N, 2) # hold shape parameters for beta distribution for
# each observation (Y_n, X_n)

## generate data
for (n in 1:N) {
    # X_n = (X_n1, X_n2, ... , X_np) ~ Unif [0.05, 0.95]^p
    X[n,] = runif(D, 0.05, 0.95) # could generate all covariates at once, 
    # but this seems a little more intuitive
    
    # store shape parameters for beta distribution
    shape_mat[n, 1] = 4 * X[n, 1] + 3 * X[n, 2]^2
    shape_mat[n, 2] = 10 * X[n, 2]
    
    # y_n ~ Beta(shape1, shape2)
    y[n] = rbeta(1, shape1 = shape_mat[n,1], shape2 = shape_mat[n, 2])
} 

# test functions (covariate-dependent) -----------------------------------------

source("initVarParams.R")

theta = initVarParams(y, X, N, D, K, max_iter)


