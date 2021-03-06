
## initPriors.R -- covariate-DEPENDENT case w/ variable selection

# create a prior list using the priors passed in for each parameter
# input:
        # y         : (N x 1) -- response values
        # X         : (N x D) -- covariates stored row-wise
        # K         : # of mixture componenets
        # m_0       : prior mean for beta_d; k = 1, ... , D
        # xi_0      : scale for precision matrix for beta_d; d = 1, ... , D
        # a_0       : prior shape parameter for tau_k; k = 1, ... , K
        # b_0       : prior rate parameter for tau_k; k = 1, ... , K
        # g_0       : prior mean for gamma_k; k = 1, ... , K
        # Sigma_0   : prior covariance for gamma_k; k = 1, ... , K
        # max_iter  : max # of iters to run CAVI
        # tol       : tolerlance for evaluating change in elbo
        # VERBOSE   : if TRUE, then detailed convergence output will be produced
# output: 
        # prior     : list containing the prior parameters

initPriors = function(y, X, K, 
                      m_0, xi_0, pi_d,             # beta components
                      a_0, b_0,                    # tau components
                      g_0, Sigma_0,                # gamma components 
                      max_iter, tol, VERBOSE) {
    
    prior = list()
    
    # response, design matrix
    prior$y        = y
    prior$X        = X
    
    # square of response
    prior$y_sq = y^2
    
    # num observations, dim of covariates, num of clusters
    prior$N        = nrow(X)
    prior$D        = ncol(X)
    prior$K        = K
    
    # prior mean and precision for beta_1:K
    # prior$m_0      = m_0
    # prior$Lambda_0 = Lambda_0
    
    # prior$Lambda0_m0    = Lambda_0 %*% m_0
    # prior$m0_Lambda0_m0 = c(quadMult(m_0,  Lambda_0))
    
    # prior mean and precision (scale) for slab component of beta_1:D 
    prior$m_0      = m_0
    prior$xi_0     = xi_0
    
    # probability of beta_d drawn from the slab (choose all to be same prob)
    prior$pi_d     = rep(pi_d, D)
    prior$log_odds = logit(prior$pi_d)
    
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















# end of initPriors.R 
