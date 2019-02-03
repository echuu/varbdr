
## set to cov_dep working directory, functions are named identically
## in the covariate independent directory
setwd("C:/Users/chuu/varbdr/code/cov_dep")


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

theta$alpha
theta$xi
theta$lambda

theta$beta_k
theta$tau_k
theta$gamma_k

theta$a_k
theta$b_k

n = 30
params = list(shape1 = shape_mat[n,1], shape2 = shape_mat[n,2])
p1 = plotDensities(y_grid, X[n,], dbeta, params, p_y, theta, prior, K)
p1$overlay

# compare iterative changes for the density curve

# true density
beta_n = stat_function(aes(x = y_grid, y = ..y..), 
                       fun = dbeta, colour = 'black', n = 500, size = 1,
                       args = list(shape1 = shape_mat[n,1], 
                                   shape2 = shape_mat[n,2])) 

curve_2 = data.frame(x = y_grid, y = theta$dc[[2]][n,])
py_plot2 = geom_line(aes(x = curve_2$x, y = curve_2$y), 
                     colour = "red", size = 0.9)

curve_50 = data.frame(x = y_grid, y = theta$dc[[3]][n,])
py_plot50 = geom_line(aes(x = curve_50$x, y = curve_50$y), 
                      colour = "orange", size = 0.9)


curve_100 = data.frame(x = y_grid, y = theta$dc[[4]][n,])
py_plot100 = geom_line(aes(x = curve_100$x, y = curve_100$y), 
                       colour = "yellow", size = 0.9)

curve_150 = data.frame(x = y_grid, y = theta$dc[[5]][n,])
py_plot150 = geom_line(aes(x = curve_150$x, y = curve_150$y), 
                      colour = "green", size = 0.9)

curve_151 = data.frame(x = y_grid, y = theta$dc[[6]][n,])
py_plot151 = geom_line(aes(x = curve_151$x, y = curve_151$y), 
                       colour = "blue", size = 0.9)

curve_152 = data.frame(x = y_grid, y = theta$dc[[7]][n,])
py_plot152 = geom_line(aes(x = curve_152$x, y = curve_152$y), 
                       colour = "purple", size = 0.9)

p + py_plot2 + py_plot50 + py_plot100 + py_plot150 + py_plot151 + 
    py_plot152 + beta_n


