


## mStep.R -- covariate dependent with sparsity assumption
## perform one iteration of the variational m-step



# input: 
#          theta : list of variational parameters
# output:
#          theta : list of variational parameters with 
#                  variational parameters updated

mStep = function(prior, theta) {
    
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
    # theta = qgamma$theta
    
    
    # (3) update q(beta, omega) = q(beta_1, omega_1) x ... x q(beta_D, omega_D)
    qbeta = spikeSlabUpdate(prior, theta)
    
    # TODO: stick these updates in spikeSlabUpdate() function (see (2))
    theta$Q_d      = qbeta$Q_d
    theta$Q_d_inv  = qbeta$Q_d_inv      # V [ beta_d | omega_d = 1 ] = Q_d^(-1)
    theta$m_d      = qbeta$m_d          # E [ beta_d | omega_d = 1 ] = m_d
    theta$lambda_d = qbeta$lambda_d     # posterior inclusion probabilities
    
    
    # (3.1) update E[beta_k] (D x 1), V(beta_k) (D x D)
    # col-wise expectation, covariance for covariates (used in tau_k updates)
    # each updated with the posterior means, variances under q(beta_d)
    
    # TODO: stick this into the outer loop of Sigma_k construction after tested
    for (k in 1:K) {
        theta$m_k[,k] = theta$lambda_d * theta$m_d[k,] # (D x 1) column vector
    }
    
    # construct Sigma_k (D x D) for each of the beta_k's
    for (k in 1:K) {
        
        for (d in 1:D) { # compute the D diagonal entries of Sigma_k
            l_d = theta$lambda_d[d]
            md_md_t = tcrossprod(theta$m_d[,d], theta$m_d[,d])
            
            theta$Sigma_k[d, d, k] = l_d * (theta$Q_d_inv[,,d] + 
                                                (1 - l_d) * md_md_t)[k,k]  
        } # end inner for()
        
    } # end outer for()
    
    
    # TODO: compute these outside of the m-step in the main CAVI loop
    # update model parameters w/ posterior means: ------------------------------
    
    # unconditional posterior mean, E[beta_d]
    theta$beta_d = theta$lambda_d * theta$m_d
    
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

