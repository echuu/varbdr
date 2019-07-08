


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
    
    # variational distributions to be updated: q(tau), q(gamma), q(beta, omega)
    
    # (1) update: q(tau) = q(tau_1) x ... x q(tau_K)
    theta = precisionUpdate(prior, theta)
    
    # (2) update: q(gamma) = q(gamma_1) x ... x q(gamma_K)
    theta = weightUpdate(prior, theta)
    
    # (3) update q(beta, omega) = q(beta_1, omega_1) x ... x q(beta_D, omega_D)
    theta = spikeSlabUpdate(prior, theta)
    
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
    
    
    return(theta)
    
}

