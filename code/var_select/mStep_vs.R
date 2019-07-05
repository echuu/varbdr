


## mStep.R -- covariate dependent with sparsity assumption
## perform one iteration of the variational m-step



# input: 
#          theta : list of variational parameters
# output:
#          theta : list of variational parameters with 
#                  variational parameters updated

mStep = function(theta_vs, prior_vs) {
    
    # extract design, response, and dimensions for local use
    X = prior$X                                 # (N x D) design matrix
    y = prior$y
    N = prior$N                                 # num of observations
    K = prior$K                                 # num of clusters
    D = prior$D                                 # dimension of covariates
    
    # commonly used matrices
    I_D        = diag(1, D)                      # (D x D) : identity matrix
    I_K        = diag(1, K)                      # (K x K) : identity matrix
    X_mu       = X %*% theta$mu_k                # (N x K) : (N x D) * (D x K)
    
    
    # variational distributions to be updated:
    # q(tau), q(gamma), q(beta, omega)
    # (1) update: q(tau) = q(tau_1) x ... x q(tau_K)
    # (2) update: q(gamma) = q(gamma_1) x ... x q(gamma_K)
    
    # (3) update q(beta, omega) = q(beta_1, omega_1) x ... x q(beta_D, omega_D)
    
    
    
    
    
    
    
    
    
    
    
    
    
}

