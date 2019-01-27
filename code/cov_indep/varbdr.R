

# varbdr.R -- covariate-INDEPENDENT case

source("initPriors.R")
source("initVarParams.R")
source("eStep.R")
source("mStep.R")
source("elbo.R")
source("misc.R")

## conditional density estimation using mixture of experts with covariate 
## INDEPENDENT weights + VB for faster inference

## initialize model parameters -------------------------------------------------


## globals



## input:
#         y        : (N x 1) -- response values
#         X        : (N x D) -- covariates stored row-wise
#         K        : # of mixture componenets
#         alpha_0  : prior concentration parameter for pi_k
#         m0       : prior mean for beta_k; k = 1,...,K
#         Lambda0  : prior (scaled) precision matrix for beta_k; k = 1,...,K
#         a0       : prior shape parameter for tau_k; k = 1,...,K
#         b0       : prior rate parameter for tau_k; k = 1,...,K
varbdr = function(y, X, K = 3, 
                  alpha_0 = rep(1 / K, K),                      # dir param
                  m_0 = c(colMeans(X)),                         # normal params
                  Lambda_0 = diag(rep(1, ncol(X))), 
                  a_0 = 1, b_0 = 1,                             # gamma params
                  max_iter = 15, tol = 1e-4, VERBOSE = TRUE) {
    
    # TODO: set default values for a_0, b_0 
    
    
    X = as.matrix(X)         # N X D design matrix of covariates
    D = NCOL(X)              # Number of features
    N = NROW(X)              # Number of observations
    L = rep(-Inf, max_iter)  # Store the variational lower bounds
    
    # prior specification: 
    #           pi ~ Dir (alpha_0)
    # mu_k | tau_k ~ N  (m0, (tau_k * Lambda0)^{-1})
    #        tau_k ~ Ga (a0, b0)
    
    # intialize prior object using the prior parameters passed in
    prior = initPriors(y, X, K, alpha_0, m_0, Lambda_0, a_0, b_0, 
                       max_iter, tol, VERBOSE)
    
    # initialize variational parameters
    #     N = number of observations
    #     D = dimension of covariaftes
    #     K = number of clusters
    theta = initVarParams(y, X, N, D, K, max_iter)
    
    # begin CAVI ---------------------------------------------------------------
    
    # perform coordinate ascent on variational parameters until EITHER:
    #     (1) ELBO converges
    #     (2) max_iter CAVI iterations have been performed
    for (i in 2:max_iter) {
        
        # e-step, m-step
        theta = eStep(theta, prior)
        theta = mStep(theta, prior)
        
        # compute ELBO
        theta$L[i] = elbo(theta, prior)
        
        ## some calculation is done here in cov-independent code to recompute
        ## 2 expectations
        
        # check for convergence
        if (checkELBO(theta, prior)) {    # checkElbo() in misc.R
            break
        }
        
    } # end of CAVI loop -------------------------------------------------------
    
    # return list containing udpated variational parameters
    
    return(theta)
    
} # end of varbdr() function


