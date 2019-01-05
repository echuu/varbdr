

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
varbdr = function(X, K = 3, m_0 = c(colMeans(X)), Lambda0 = I_D, 
                  a0, b0, g0 = 0, Sigma0 = I_D, max_iter = 500, 
                  tol = 1e-4, is_animation = FALSE, VERBOSE = FALSE) {
    
    X = as.matrix(X)         # N X D design matrix of covariates
    D = NCOL(X)              # Number of features
    N = NROW(X)              # Number of observations
    L = rep(-Inf, max_iter)  # Store the variational lower bounds
    
    # prior specification: 
        # mu_k | tau_k ~ N  (m0, (tau_k * Lambda0)^{-1})
        #        tau_k ~ Ga (a0, b0)
        # gamma_k      ~ N  (g0, Sigma0)
    
    # store the priors into a prior list
    prior = list()

    # initialize variational parameters    
    theta = vb_initVarParams(X, K)

    
    # in gmVB_0.R code, we compute two of the three expectations outside of the
    # cavi loop since these are normally computed at the END of each cavi
    # iteration, so we need these to have values for the FIRST iteration
    
    
    
    # begin CAVI ---------------------------------------------------------------
    
    # perform coordinate ascent on variational parameters until EITHER:
        # (1) ELBO converges
        # (2) max_iter iterations have been performed
    for (i in 2:max_iter) {
        
        # Variational E-step ---------------------------------------------------
        
        # update responsibilities r_nk; involves E[ln tau_k], E[gamma_k],
        # E[tau_k * (y_n - x_n'beta_k)^2], E[ln sum_j (x_n'gamma_j)]
        
        for (k in 1:K) {
            
            
            
        } # end of updating rho_nk
        
        # update E[z_nk] = r_nk ------------------------------------------------
        # log of the normalizing constant for the rho_nk's
        # Z = log { sum_{j=1}^{k} exp( ln rho_{nj} ) }
        logZ     = apply(log_rho_nk, 1, log_sum_exp) # log of norm. constant
        log_r_nk = log_rho_nk - logZ                 # log of r_nk
        r_nk     = apply(log_r_nk, 2, exp)           # exponentiate each column
        
        # this is done within the estep() function
    
    
    
        
        # Variational M-step ---------------------------------------------------
        
        
    } # end of CAVI
    
    
    
    # prepare output -----------------------------------------------------------
    
    # return list containing udpated variational parameters
    
    out = list()
    
    return(out)
    
} # end of varbdr() function


