

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
    
    # prior specification: 
        # mu_k | tau_k ~ N  (m0, (tau_k * Lambda0)^{-1})
        #        tau_k ~ Ga (a0, b0)
        # gamma_k      ~ N  (g0, Sigma0)
    
    
    
    # initialize storage for variational parameters ----------------------------
    
    X = as.matrix(X)         # N X D design matrix of covariates
    D = NCOL(X)              # Number of features
    N = NROW(X)              # Number of observations
    L = rep(-Inf, max_iter)  # Store the variational lower bounds
    
    
    r_nk       = matrix(0, N, K) # N x K : normalized responsibilities
                                 #         rho_nk / sum_j rho_nj
                                 # can calculate this by exponentiating 
                                 # log_rho_nk and dividing each element by the
                                 # sum of the ROW it belongs in
    log_r_nk   = matrix(0, N, K) # N x K : log(r_nk)
    log_rho_nk = matrix(0, N, K) # N x K : log(rho_nk)
    
    
    # variational parameters used to compute Q_k_inv, eta_k via Bouchard
    alpha   = numeric(N)            # N x 1 : alpha_n used to compute xi_{n,1:K}
    xi      = matrix(0, N, K)       # N x K : nth row used to compute alpha_n
    
    # variational parameters for gamma_1:K ~ N(gamma_k | mu_k, Q_k^{-1})
    Q_k_inv = array(0, c(D, D, K))   # K x (D x D) precisions for gamma_k
    eta_k   = matrix(0, D, K)        # D x K       Q_k_inv * eta_k = mu_k
    mu_k    = matrix(0, D, K)        # D x K       mean of gamma_k
    
    # variational parameters for beta_k | tau_k
    V_k_inv = array(0, c(D, D, K))   # K x (D x D) scaled precision for beta_k
    zeta_k  = matrix(0, D, K)        # D x K       V_k_inv * zeta_k = m_k
    m_k     = matrix(0, D, K)        # D x K       mean of gamma_k
    
    
    # precompute values that are used regularly in CAVI ------------------------
    # values that are used in the r_nk update
    
    # begin CAVI ---------------------------------------------------------------
    
    # perform coordinate ascent on variational parameters until EITHER:
        # (1) ELBO converges
        # (2) max_iter iterations have been performed
    for (i in 2:max_iter) {
        
    
        # Variational E-step ---------------------------------------------------
        
    
    
    
    
        # Variational M-step ---------------------------------------------------
        
        
    } # end of CAVI
    
    
    
    # prepare output -----------------------------------------------------------
    
    # return list containing udpated variational parameters
    
    out = list()
    
    return(out)
    
} # end of varbdr() function


