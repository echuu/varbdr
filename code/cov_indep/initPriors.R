
## initPriors.R -- covariate-INDEPENDENT case

# create a prior list using the priors passed in for each parameter
# input:
#         y        : (N x 1) -- response values
#         X        : (N x D) -- covariates stored row-wise
#         K        : # of mixture componenets
#         alpha_0  : prior concentration parameter for pi; k = 1, ... , K
#         m_0      : prior mean for gamma_k; k = 1, ... , K
#         Lambda0  : prior (scaled) precision matrix for beta_k; k = 1, ... , K
#         a_0      : prior shape parameter for tau_k; k = 1, ... , K
#         b_0      : prior rate parameter for tau_k; k = 1, ... , K
#         max_iter  : max # of iters to run CAVI
#         tol       : tolerlance for evaluating change in elbo
#         VERBOSE   : if TRUE, detailed convergence output will be produced
# output: 
#         prior    : list containing the prior parameters
initPriors = function(y, X, K, alpha_0, m_0, Lambda_0, a_0, b_0, 
                      max_iter, tol, VERBOSE) {
    
    prior = list()
    
    # response vector (N x 1), design matrix (N x K)
    prior$y = y    # (N x 1)
    prior$X = X    # (N x K)
    
    # num of observations, dimension of covariates, num of clusters
    prior$N = nrow(X)
    prior$D = ncol(X)
    prior$K = K
    
    
    # prior concentration parameter for alpha_1, ... , alpha_K
    prior$alpha_0 = alpha_0 # (K x 1)
    
    # prior mean and (scaled) precision matrix for beta_1, ... , beta_K
    prior$m_0      = m_0         # (D x 1)
    prior$Lambda_0 = Lambda_0    # (D x D)
    
    # prior shape, rate for tau_1, ... ,tau_K
    prior$a_0 = a_0    # (1 x 1)
    prior$b_0 = b_0    # (1 x 1)
    
    # other algorithm-related parameters that are determined before CAVI iters
    prior$max_iter = max_iter
    prior$tol      = tol
    prior$VERBOSE  = VERBOSE
    
    return(prior)
    
} # end initPriors() function

