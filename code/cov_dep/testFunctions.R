
## testFunctions.R -- test/troubleshoot functions (covariate-DEPENDENT)


# generate data ----------------------------------------------------------------

N        = 100  # number of observations
D        = 2    # number of covariates 
K        = 5    # number of clusters
max_iter = 200  # max number of iterations 

tol = 1e-4
VERBOSE = TRUE


# generate X1,...,Xp ~ Unif[0.05, 0.95]
## initialize storage for response, covariates
y = numeric(N)       # (N x 1) vector of response variables
X = matrix(0, N, D)  # (N x D) matrix; covariates for y_n are stored row-wise  

shape_mat = matrix(0, N, 2) # hold shape parameters for beta distribution for
# each observation (Y_n, X_n)

## generate data ---------------------------------------------------------------
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
source("initPriors.R")
source("initVarParams.R")

## prior parameters
m_0 = c(colMeans(X))                                    
Lambda_0 = diag(rep(1, ncol(X)))
a_0 = 1
b_0 = 1                             
g_0 = rep(0, ncol(X))
Sigma_0 = diag(rep(1, ncol(X)))


# initPriors() function tested and works
prior = initPriors(y, X, K, m_0, Lambda_0, a_0, b_0, g_0, Sigma_0,
                   max_iter, tol, VERBOSE)

# initVarParams() function tested and works
theta = initVarParams(y, X, N, D, K, max_iter)


source("eStep.R")
theta_estep = eStep(theta, prior)

source("mStep.R")
theta = mStep(theta, prior)

# currently testing functions below --------------------------------------------

source("elbo.R")
theta$L[i] = elbo(theta, prior)



X_mu = X %*% theta$mu_k
t(X_mu[n,]) %*% theta$lambda[n,]


r_vec = rnorm(N)
l_vec = rnorm(N)
a_vec = rnorm(N)
    
mat_result = t(X) %*% (r_vec * (0.5 + 2 * (l_vec * a_vec)))

loop_res = 0
for (n in 1:N) {
    loop_res = loop_res + r_vec[n] * (0.5 + 2 * l_vec[n] * a_vec[n]) * X[n,]
}

loop_result = 0
for (k in 1:K) {
    loop_result = loop_result + t(theta$m_k[,k]) %*% theta$m_k[,k]
}


x_grid = seq(0.1, 100, length = 1000)
lambda_x = lambda_xi(x_grid)
xl_df = data.frame(x = x_grid, y = lambda_x)
ggplot(xl_df, aes(x, y)) + geom_point(size = 0.7)




