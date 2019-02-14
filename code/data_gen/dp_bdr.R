
# dp_bdr.R -- contains the data generating functions in Dunson and Park paper:
#             1-d examples



# r_dpmix1()     : generate data according to the first simulation example in
#                  Bayesian Density Regression (Dunson, Park)
# input:  
#          N     : number of observations/data points to generate
#          seed  : seed for RNG 
# output: 
#          y     : (N x 1) vector of responses
#          X     : (N x 2) design matrix, covariates stored row-wise
#          shape : (N x 2) shape parameters for each y_n
r_dpmix1 = function(N, seed = 100) {
    set.seed(seed)
    sigma_sq = 0.01
    x_i1 = rep(1, N) # 1st column of X
    x_i2 = runif(N)  # 2nd column of X
    
    X = matrix(c(x_i1, x_i2), nrow = N, byrow = FALSE)
    mu_y = -1 + 2 * X[,2] # (N x 1) vector of conditional means
    
    y = rnorm(N, mean = mu_y, sd = sqrt(sigma_sq)) # (N x 1) vector of responses
    
    return(list(y = y, X = X, param_vec = mu_y))
}


# r_dpmix2()     : generate data according to the second simulation example in
#                  Bayesian Density Regression (Dunson, Park)
# input:  
#          N     : number of observations/data points to generate
#          seed  : seed for RNG 
# output: 
#          X     : (N x 2) design matrix, covariates stored row-wise
#          y     : (N x 1) vector of responses
#          shape : (N x 2) shape parameters for each y_n
r_dpmix2 = function(N, seed = 100) {
    
    set.seed(seed)
    
    x_i1 = rep(1, N) # 1st column of X
    x_i2 = runif(N)  # 2nd column of X
    
    X = matrix(c(x_i1, x_i2), nrow = N, byrow = FALSE)
    
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
    
    
    return(list(y = y, X = X, param_mat = data.frame(mu_1, mu_2)))
}




# end dp_bdr.R
