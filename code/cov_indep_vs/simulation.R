
# simulation.R --- generate (x,y) data for bayesian density regression

# gmm_sim_xy() --- generate data points (x, y) where y is drawn from a 
#                  conditional density f(y|x)
# input: 
#         N    :  # of observations to generate
#         D    :  dimension of the covariate
# output: list containing X, y 
#         X    :  (N x D) design matrix
#         y    :  (N x 1) response vector
gmm_sim_xy = function(N, D, seed = 100) {
    
    set.seed(seed)
    
    
    X = matrix(runif(N * D), N, D)
    
    # means of the two gaussian components --> X[,2:D] are noise covariates
    # which are ideally zero'd out when implementing the sparsity assumption
    mu_1 = X[,1]
    mu_2 = -2 * X[,1]
    
    # variances of the two gaussian components
    sigsq_1 = 0.01
    sigsq_2 = 0.04
    
    # probability of being drawn from each mixture
    p = 0.3                     # prob of drawing from N ( y | mu_1, sigsq_1)
    
    # 0 if draw from N ( y | mu_1, sigsq_1 ), 1 if from N ( y | mu_2, sigsq_2 ) 
    mix_indic = (runif(N) >= p) 
    
    N_mix1 = sum(!mix_indic)   # num of draws from N ( y | mu_1, sigsq_1 )
    N_mix2 = sum(mix_indic)    # num of draws from N ( y | mu_2, sigsq_2 )
    
    # draw y ~ f( y | x )
    y = numeric(N)
    y[!mix_indic] = rnorm(N_mix1, mean = mu_1[!mix_indic], sd = sqrt(sigsq_1))
    y[mix_indic]  = rnorm(N_mix2, mean = mu_2[mix_indic], sd = sqrt(sigsq_2))
    
    return(list(X = X, y = y))
    
} # end of gmm_sim_xy() function
