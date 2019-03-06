

# r_binorm()     : generate data according to the bimodal conditional density
#                  example in Shape-Constrained Univariate Density Estimation,
#                  (Dasgupta, Pati, Jermyn, Srivastava)
#                  
# input:  
#          N     : number of observations/data points to generate
#          seed  : seed for RNG 
# output: 
#          X     : (N x 2) design matrix, covariates stored row-wise
#          y     : (N x 1) vector of responses
#          shape : (N x 2) shape parameters for each y_n
r_binorm = function(N, seed = 100) {
    set.seed(seed)
    
    x_vec = rnorm(N)  # 2nd column of X
    
    X = matrix(x_vec, nrow = N)
    
    # mixture parameters
    mu_1    = X - 1.5
    mu_2    = X + 1.5
    sigma_1 = 0.5
    sigma_2 = 0.5
    
    # p = exp(-0.5 * x_vec)      # probability that we draw from mixture 2
    p = 0.5
    mix_indic = runif(N) >= p  # determine which density to draw from
    
    N_mix1 = sum(!mix_indic)   # number of draws from first mixture density
    N_mix2 = sum(mix_indic)    # number of draws from second mixture density
    
    y = numeric(N)
    y[!mix_indic] = rnorm(N_mix1, mean = mu_1[!mix_indic], sd = sigma_1)
    y[mix_indic]  = rnorm(N_mix2, mean = mu_2[mix_indic],  sd = sigma_2)
    
    
    return(list(y = y, X = X, param_mat = data.frame(mu_1, mu_2)))
}




