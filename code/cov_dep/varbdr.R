

# varbdr.R

## conditional density estimation using mixture of experts with covariate 
## dependent weights + VB for faster inference


## global params ---------------------------------------------------------------

D     = 2           # dimension of covariate, 2 <= d <= 100
K     = 10          # number of clusters
N     = 1           # number of samples
I_D   = diag(1, D)  # D X D  identity matrix
df    = 100         # degrees of freedom for Wishart distribution
scale = 1           # scale for the covariance matrix


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
varbdr = function(y, X, K = 3, m_0 = c(colMeans(X)), Lambda_0 = I_D, 
                  a_0, b_0, g_0 = 0, Sigma_0 = I_D, max_iter = 500, 
                  tol = 1e-4, is_animation = FALSE, VERBOSE = TRUE) {
    
    X = as.matrix(X)         # N X D design matrix of covariates
    D = NCOL(X)              # Number of features
    N = NROW(X)              # Number of observations
    L = rep(-Inf, max_iter)  # Store the variational lower bounds
    
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
    theta = initVarParams(N, D, K)
    
    # begin CAVI ---------------------------------------------------------------
    
    # perform coordinate ascent on variational parameters until EITHER:
        # (1) ELBO converges
        # (2) max_iter iterations have been performed
    for (i in 2:max_iter) {
        
        # e-step, m-step
        theta = eStep(theta)
        theta = mStep(theta, prior)
        
        # compute ELBO
        theta$L[i] = elbo(theta, prior)
        
        ## some calculation is done here in cov-independent code to recompute
        ## 2 expectations
        
        # check for convergence
        if (elboConverged(theta, prior)) {
            break
        }
        
    } # end of CAVI loop
    
    # return list containing udpated variational parameters
    
    return(theta)
    
} # end of varbdr() function


