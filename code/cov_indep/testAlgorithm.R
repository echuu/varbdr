
## testAlgorithm.R -- run vb algorithm for covariate-INDEPENDENT case

library(ggplot2)
library(reshape2)


N = 100  # number of observations
D = 2    # number of covariates 
K = 5    # number of clusters

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

# plot each of the conditional densities (N densities)
# overlay each of the density plots
y_grid = seq(0, 1, len = 500)
p = qplot(y_grid, geom = "blank")
for (n in 1:N) {
    beta_n = stat_function(aes(x = y_grid, y = ..y..), 
                           fun = dbeta, colour = 'grey', n = 100, 
                           args = list(shape1 = shape_mat[n,1], 
                                       shape2 = shape_mat[n,2]))
    p = p + beta_n
}
p + labs(x = "y", y = "p (y | x)", 
         title = "Conditional Densities for varying values of x")

## begin VB algorithm for conditional density estimation -----------------------

run = function() {
    source("varbdr.R")              # load the CAVI algorithm for BDR
    theta = varbdr(y = y, X = X, K)    # run algorithm
    return(theta)
}


theta = run()                   # run algorithm

# variational parameters after CAVI finishes
theta$alpha_k

theta$m_k
theta$tau_k

theta$a_k
theta$b_k


source("approxDensity.R")

data = list(y = y[1], x = X[1,])
p_y(theta, prior, K, data) # approximate density

dbeta(y[1], shape_mat[1,1], shape_mat[1,2]) # true density

# generate density plots of true density overlayed with approximate density
n = 83
data_ygrid = list(y = y_grid, x = X[n,])     # covariates for the n-th observ.
p_ygrid = p_y(theta, prior, K, data_ygrid)

apx_beta_n = data.frame(x = y_grid, y = p_ygrid)

# true density plot (indexed by n in the shape parameters)
beta_n = stat_function(aes(x = y_grid, y = ..y..), 
                       fun = dbeta, colour = 'red', n = 500, size = 1,
                       args = list(shape1 = shape_mat[n,1], 
                                   shape2 = shape_mat[n,2]))  


p = ggplot(apx_beta_n, aes(x, y)) + geom_point(colour = 'blue', size = 0.9)
p + beta_n




# evalulate performance (?); this part still a little unsure of what to do
# do i need to have a predictive density ready? doesn't really make sense to 
# look at coefficients because we introduce those artificially as part of
# the model.


