
## initPriors.R

# create a prior list using the priors passed in for each parameter
# input:
        # y        : (N x 1) -- response values
        # X        : (N x D) -- covariates stored row-wise
        # K        : # of mixture componenets
        # m0       : prior mean for beta_k; k = 1,...,K
        # Lambda0  : prior (scaled) precision matrix for beta_k; k = 1,...,K
        # a0       : prior shape parameter for tau_k; k = 1,...,K
        # b0       : prior rate parameter for tau_k; k = 1,...,K
        # g0       : prior mean for gamma_k; k = 1,...,K
        # Sigma0   : prior covariance for gamma_k; k = 1,...,K
# output: 
        # prior    : list containing the prior parameters
initPriors = function(y, X, m_0, Lambda0, a_0, b_0, g0, Sigma0,
                      max_iter, tol, VERBOSE) {
    
    prior = list()
    
    # response, design matrix
    prior$y        = y
    prior$X        = X
    
    # prior mean and precision for beta_1:K
    prior$m_0      = m_0
    prior$Lambda_0 = Lambda_0
    
    # prior shape, rate for tau_1:K
    prior$a_0      = a_0
    prior$b_0      = b_0
    
    # prior mean, precision for gamma_1:K
    prior$g_0      = g_0
    prior$Sigma_0  = Sigma_0
    
    # other algorithm-related parameters, 
    prior$max_iter = max_iter
    prior$tol      = tol
    prior$VERBOSE  = VERBOSE
    
    
    return(prior)
} # end initPriors() function

# end initPriors.R
