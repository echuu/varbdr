
## testAlgorithm.R -- run vb algorithm for covariate-INDEPENDENT case

library(ggplot2)

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
set.seed(100)
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


## evaluate results ------------------------------------------------------------

# examine observation n = 32, per every 2 iteration of CAVI
n = 32
y_grid = seq(0, 1, len = 500)
iters = seq(2, theta$curr, 2)

# true density
beta_n = stat_function(aes(x = y_grid, y = ..y..), 
                       fun = dbeta, colour = 'black', n = 500, size = 1,
                       args = list(shape1 = shape_mat[n,1], 
                                   shape2 = shape_mat[n,2])) 
i = 2
i = 4
i = 6
i = 10
i = 20
i = 52

for (i in 1:length(iters)) {
    # look in the n-th row of each data frame: theta$dc
    curve_1 = data.frame(x = y_grid, y = theta$dc[[iters[1]]][n,])
    py_plot1 = geom_line(aes(x = curve_1$x, y = curve_1$y), 
                        colour = "blue", size = 0.9)
    p + py_plot + beta_n
    
    curve_26 = data.frame(x = y_grid, y = theta$dc[[iters[26]]][n,])
    py_plot26 = geom_line(aes(x = curve_26$x, y = curve_26$y), 
                        colour = "red", size = 0.9)
    p + py_plot1 + py_plot26 + beta_n
    
    
}





### testing --------------------------------------------------------------------


# look at the density plots for observation n = 63 per iteration of cavi
source("approxDensity.R")

p = qplot(y_grid, geom = "blank")

curves_approx = densityCurve(p_y, theta, prior, X, K,
                             y_grid = seq(0, 1, len = 500))

n = 33
params = list(shape1 = shape_mat[n,1], shape2 = shape_mat[n,2])
p1 = plotDensities(y_grid, X[n,], dbeta, params, p_y, theta, prior, K)
p1$p


beta_n = stat_function(aes(x = y_grid, y = ..y..), 
                       fun = dbeta, colour = 'blue', n = 500, size = 1,
                       args = list(shape1 = shape_mat[n,1], 
                                   shape2 = shape_mat[n,2])) 

for (n in 1:N) {
    py_df = data.frame(x = y_grid, y = curves_approx[n,])
    
    p + geom_line(aes(x = py_df$x, y = py_df$y), colour = "red", size = 0.9) +
        beta_n
}

# ------------------------------------------------------------------------------


# generate density plots of true density overlayed with approximate density

source("approxDensity.R")


## automation of density evaluations, overlay densities way --------------------
n = 4
params = list(shape1 = shape_mat[n,1], shape2 = shape_mat[n,2])
p1 = plotDensities(y_grid, X[n,], dbeta, params, p_y, theta, prior, K)
p1


# end of testAlgorithm.R file
