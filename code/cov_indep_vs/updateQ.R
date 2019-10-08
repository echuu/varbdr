
# updateQ.R --- implementation of the variational updates for covariate
#               independent weights with sparsity assumption in the gaussian
#               components


# updateQ(): performs one update for each of the variational distributions
#   (0) q(Z)
#   (1) q(pi)
#   (2) q(tau)
#   (3) q(beta, omega)
updateQ = function(prior, theta) {
    
    # (0) update q(Z)
    theta = rnkUpdate(prior, theta)
    
    # (1) update q(pi) = Dir(pi | alpha_1, ... , alpha_K)
    theta = wtUpdate(prior, theta)
    
    # (2) update q(tau) = Ga(tau_1 | a_1, b_2) x ... x Ga(tau_K | a_K, b_K)
    theta = precUpdate(prior, theta)
    
    # (3) update q(beta, omega)
    # (3.1) update q(beta_d | omega_d = 1) = N(beta_d | m_d, Q_d^{-1}), d = 1:D
    # (3.2) update omega_d, d = 1:D
    theta = ssUpdate(prior, theta)
    
    # (3.3) update the unconditional mean and variance of beta_k
    theta = betak_update(prior, theta)
    
    # update variables used in the final model
    theta = updateModelParams(prior, theta)
    
    return(theta)
    
} # end updateQ() function

# end upadateQ.R
