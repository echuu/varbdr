
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
r_f2 = function(N, X, seed = 100) {
    
    set.seed(seed)
    
    # mixture parameters
    mu_1    = X[,2]
    mu_2    = X[,2]^4
    sigsq_1 = 0.01
    sigsq_2 = 0.04
    
    p = exp(-2 * X[,2])        # probability of drawing from mixture 1
    mix_indic = runif(N) >= p  # determine which density to draw from
    
    N_mix1 = sum(!mix_indic)   # number of draws from first mixture density
    N_mix2 = sum(mix_indic)    # number of draws from second mixture density
    
    y = numeric(N)
    y[!mix_indic] = rnorm(N_mix1, mean = mu_1[!mix_indic], sd = sqrt(sigsq_1))
    y[mix_indic]  = rnorm(N_mix2, mean = mu_2[mix_indic],  sd = sqrt(sigsq_2))
    
    
    return(list(y = y, mix = mix_indic))
}

# give conditional density of y | x
d_rf2 = function(x, N = 1000) {
    
    mu_1    = x
    mu_2    = x^4
    sigsq_1 = 0.01
    sigsq_2 = 0.04
    
    N_evals = length(y_grid)
    
    # p = exp(-2 * X)                  # probability of drawing from mixture 1
    #mix_indic = runif(N_evals) >= p  # determine which density to draw from
    
    N_mix1 = sum(!mix_indic)         # num of draws from first mixture density
    N_mix2 = sum(mix_indic)          # num of draws from second mixture density
    
    y = numeric(N)
    y[!mix_indic] = rnorm(N_mix1, mean = mu_1, sd = sqrt(sigsq_1))
    y[mix_indic]  = rnorm(N_mix2, mean = mu_2, sd = sqrt(sigsq_2))

    
    f_y = numeric(N_evals)
    f_y[!mix_indic] = dnorm(y[!mix_indic], mean = mu_1, sd = sqrt(sigsq_1))
    f_y[mix_indic]  = dnorm(y[mix_indic],  mean = mu_2, sd = sqrt(sigsq_2))
    
    return(data.frame(y, f_y))
}


d_rf2 = function(x, y_grid) {
    
    mu_1    = x
    mu_2    = x^4
    sigsq_1 = 0.01
    sigsq_2 = 0.04
    
    N_evals = length(y_grid)
    
    p = exp(-2 * x)                  # probability of drawing from mixture 1
    #mix_indic = runif(N_evals) >= p  # determine which density to draw from
    
    f_y = numeric(N_evals)
    f_y = p * dnorm(y_grid, mean = mu_1, sd = sqrt(sigsq_1)) + 
        (1 - p) * dnorm(y_grid, mean = mu_2, sd = sqrt(sigsq_2))
    
    return(data.frame(y = y_grid, f_y))
}


library(ggplot2)

set.seed(1)
N = 500
p = 2

# first simutlation ------------------------------------------------------------

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

df_yx = data.frame(x = X[,2], y = y)

ggplot(df_yx, aes(x, y)) + geom_point()


# second simulation ------------------------------------------------------------


y_rf2 = r_f2(N, X)
y2 = y_rf2$y
y2_mix = y_rf2$mix

# 6th sublot in Figure 3
df_yx = data.frame(x = X[,2], y = y2)
ggplot(df_yx, aes(x = y, y = fy)) + geom_point() + stat_smooth() + theme_bw()


# 1st subplot
y_grid = seq(-0.5, 1.5, length.out = 500)
x = 0.75
f_yx = d_rf2(x, y_grid)
# df_y = data.frame(x = y_grid, y = f_yx)
ggplot(f_yx, aes(x = y, y = f_y)) + geom_point(size = 0.9) 






# end dp_bdr.R
