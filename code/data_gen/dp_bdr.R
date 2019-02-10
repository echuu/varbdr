
# dp_bdr.R -- contains the data generating functions in Dunson and Park paper:
#             1-d examples

# dp() -- generate data according to the first simulation example in
#         Bayesian Density Regression (Dunson, Park)
# input:  
#          N     : number of observations/data points to generate
#          seed  : seed for RNG 
# output: 
#          X     : (N x 2) design matrix, covariates stored row-wise
#          y     : (N x 1) vector of responses
#          shape : (N x 2) shape parameters for each y_n
dp = function(N, seed = 100) {
    
    # x_i = (1, x_{i2}), x_{i2} ~ unif(0, 1)
    # f(y_i | x_i) = N(y_i; - 1 + 2x_{i2}, 0.01)
    
    
}


library(ggplot2)

set.seed(1)
N = 500
p = 2
sigma_sq = 0.01


x_i1 = rep(1, N) # 1st column of X
x_i2 = runif(N)  # 2nd column of X

X = matrix(c(x_i1, x_i2), nrow = N, byrow = FALSE)
mu_y = -1 + 2 * X[,2] # (N x 1) vector of conditional means

y = rnorm(N, mean = mu_y, sd = sqrt(sigma_sq)) # (N x 1) vector of responses


# fig 2 in DP : predictive density y_{n+1} at the 10, 25, 50, 75, 90 percentiles
#               of the empirical distribution of x_{i2}

# plot f(y|x), the true densities for varying values of x (0.1, 0.25, 0.5, 0.75)
#              to gain understanding of how the conditional density varies
#              as x goes from low to high
y_grid = seq(-1.5, 1.5, length.out = 500)

x = 0.90
f_yx = dnorm(y_grid, mean = -1 + 2 * x, sd = sqrt(sigma_sq))

df_y = data.frame(x = y_grid, y = f_yx)
ggplot(df_y, aes(x, y)) + geom_point(size = 0.9) +
    scale_x_continuous(breaks = seq(-1.5, 1.5, by = 0.5))





# end dp_bdr.R
