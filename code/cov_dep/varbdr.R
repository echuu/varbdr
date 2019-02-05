

# varbdr.R -- covariate-DEPENDENT case
## conditional density estimation using mixture of experts with covariate 
## DEPENDENT weights + VB for faster inference

HOME_DIR    = "~/varbdr/code"              # linux
HOME_DIR    = "C:/Users/chuu/varbdr/code"  # windows
COV_DEP_DIR = paste(HOME_DIR, "/cov_dep", sep = '')

setwd(COV_DEP_DIR)

source("initPriors.R")
source("initVarParams.R")
source("eStep.R")
source("mStep.R")
source("elbo.R")
source("misc.R")
source(paste(HOME_DIR, "/density.R", sep = ''))


## initialize model parameters -------------------------------------------------

## input:
    # y        : (N x 1) -- response values
    # X        : (N x D) -- covariates stored row-wise
    # K        : # of mixture componenets
    # m0       : prior mean for beta_k; k = 1,...,K
    # Lambda0  : prior (scaled) precision matrix for beta_k; k = 1,...,K
    # a0       : prior shape parameter for tau_k; k = 1,...,K
    # b0       : prior rate parameter for tau_k; k = 1,...,K
    # g0       : prior mean for gamma_k; k = 1,...,K
    # Sigma0   : prior covariance for gamma_k; k = 1,...,K
varbdr = function(y, X, K = 3, 
                  m_0 = c(colMeans(X)),                         # normal params              
                  Lambda_0 = diag(rep(1, ncol(X))), 
                  a_0 = 1, b_0 = 1,                             # gamma params
                  g_0 = 0, Sigma_0 = diag(rep(1, ncol(X))),     # normal params
                  max_iter = 1000, tol = 1e-4, VERBOSE = TRUE) {
    
    X = as.matrix(X)         # N X D design matrix of covariates
    D = NCOL(X)              # Number of features
    N = NROW(X)              # Number of observations
    # L = rep(-Inf, max_iter)  # Store the variational lower bounds
    
    # prior specification: 
        # mu_k | tau_k ~ N  (m0, (tau_k * Lambda0)^{-1})
        #        tau_k ~ Ga (a0, b0)
        # gamma_k      ~ N  (g0, Sigma0)
    
    # intialize prior object using the prior parameters passed in
    prior = initPriors(y, X, K, m_0, Lambda_0, a_0, b_0, g_0, Sigma_0,
                       max_iter, tol, VERBOSE)

    # initialize variational parameters
        # N = number of observations
        # D = dimension of covariaftes
        # K = number of clusters
    theta = initVarParams(y, X, N, D, K, max_iter)
    
    # begin CAVI ---------------------------------------------------------------
    
    # perform coordinate ascent on variational parameters until EITHER:
        # (1) ELBO converges
        # (2) max_iter iterations have been performed
    for (i in 2:max_iter) {
        
        # e-step, m-step
        theta = eStep(theta, prior)
        theta = mStep(theta, prior)
        
        # compute ELBO
        theta$L[i] = elbo(theta, prior)
        
        # obtain current density curves using each of the covariates -> N curves
        #     note: dc[[iter]] --> (N x length of grid) dataframe
        #           to get the sequential changes for the n-th iteration, 
        #           have to search along the n-th row dc[[i]], i in [2, curr]
        # theta$dc[[i]] = densityCurve(py_bouch, theta, X, K)
        
        
        # check for convergence
        if (checkELBO(theta, prior)) {
            break
        }
        
    } # end of CAVI loop
    
    # return list containing udpated variational parameters
    
    return(theta)
    
} # end of varbdr() function


