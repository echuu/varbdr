
# varbdr.R -- covariate-INDEPENDENT case

setwd("~/varbdr/code/cov_indep") # linux directory, change for windows

source("initPriors.R")
source("initVarParams.R")
source("eStep.R")
source("mStep.R")
source("elbo.R")
source("misc.R")
source("~/varbdr/code/density.R")

library(ggplot2)

## conditional density estimation using mixture of experts with covariate 
## INDEPENDENT weights + VB for faster inference

## input: 
#         y        : (N x 1) -- response values
#         X        : (N x D) -- covariates stored row-wise
#         K        : # of mixture componenets
#         alpha_0  : prior concentration parameter for pi_k
#         m0       : prior mean for beta_k; k = 1,...,K
#         Lambda0  : prior (scaled) precision matrix for beta_k; k = 1,...,K
#         a0       : prior shape parameter for tau_k; k = 1,...,K
#         b0       : prior rate parameter for tau_k; k = 1,...,K
varbdr = function(y, X, K = 4, 
                  alpha_0 = rep(1 / K, K),                      # dir param
                  m_0 = c(colMeans(X)),                         # normal params
                  Lambda_0 = diag(rep(1, ncol(X))), 
                  a_0 = 1, b_0 = 1,                             # gamma params
                  max_iter = 600, tol = 1e-4, VERBOSE = TRUE) {
    
    
    
    
    # TODO: set default values for a_0, b_0 
    X = as.matrix(X)         # N X D design matrix of covariates
    D = NCOL(X)              # Number of features
    N = NROW(X)              # Number of observations
    # L = rep(-Inf, max_iter)  # Store the variational lower bounds
    
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
    
    
    # generate sequence of y-values to be evaluated using the conditional 
    # density estimate
    # n = 30 # observation for which we want to look at the density plots
    
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
        
        # generate conditional density plot using current variational parameters
        # eventually we can just pass in all the rows of X and generate a 
        # density curve for each of the observations (maybe every 5 iterations)
        
        # obtain current density curves using each of the covariates -> N curves
        #     note: dc[[iter]] --> (N x length of grid) dataframe
        #           to get the sequential changes for the n-th iteration, 
        #           have to search along the n-th row dc[[i]], i in [2, curr]
        theta$dc[[i]] = densityCurve(p_y, theta, prior, X, K)
        
        # check for convergence
        if (checkELBO(theta, prior)) {    # checkElbo() in misc.R
            break
        }
        
    } # end of CAVI loop -------------------------------------------------------
    
    # return list containing udpated variational parameters
    
    return(theta)
    
} # end of varbdr() function


