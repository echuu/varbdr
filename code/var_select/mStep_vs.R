


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
    
    
    # variational distributions to be updated: q(tau), q(gamma), q(beta, omega)
    
    # (1) update: q(tau) = q(tau_1) x ... x q(tau_K)
    
    qtau = precisionUpdate(prior, theta)
    
    theta$a_k = qtau$a_k
    theta$b_k = qtau$b_k
    
    # (2) update: q(gamma) = q(gamma_1) x ... x q(gamma_K)
    qgamma = weightUpdate(prior, theta)
    
    theta$alpha   = qgamma$alpha
    theta$xi      = qgamma$xi
    theta$lambda  = qgamma$lambda
    theta$phi     = qgamma$phi
    theta$Q_k     = qgamma$Q_k
    theta$Q_k_inv = qgamma$Q_k_inv
    theta$mu_k    = qgamma$mu_k
    theta$eta_k   = qgamma$eta_k
    
    
    # (3) update q(beta, omega) = q(beta_1, omega_1) x ... x q(beta_D, omega_D)
    qbeta = spikeSlabUpdate(prior, theta)
    
    theta$Q_d     = qbeta$Q_d
    theta$Q_d_inv = qbeta$Q_d_inv
    theta$m_d     = qbeta$m_d
    
    
    # update model parameters w/ posterior means: ------------------------------
    
    # TODO: compute posterior mean
    theta$beta_d = theta$beta_d # to be updated
    
    # beta_k (column formulation of beta_d)
    
    # omega_d (are these used? omega_d ~ ber(lambda_d))
    
    # lambda_d (inverse logit function? something ike that)
    
    # tau_k
    theta$tau_k = theta$a_k / theta$b_k
    
    # gamma_k
    theta$gamma_k = theta$mu_k
    
    # weight not updated since it is a non-trivial calculation that really only
    # needs to be computed at the end of the algorithm to obtain the 
    # mixture density
    
    
    return(theta)
    
}

