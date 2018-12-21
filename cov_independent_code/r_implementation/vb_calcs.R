

source("misc.R")

# variational E-Step -- update r_nk using the re-estimated variational
#                       parameters from the previous iteration's m-step
# input:
#       D          : (1 x 1) dimension of covariates
#       X          : (N x D) design matrix
#       m_k        : (D x K) mean of mu_k (stored column-wise)
#       beta_k     : (K x 1) scale of precision matrix 
#       nu_k       : (K x 1) degrees of freedom for wishart distribution
#       log_pi     : (K x 1) E [ log(pi) ]
#       log_lambda : (K x 1) E [ log(det(Lambda_k)) ]
#       ** the previous two expectations are updated in the previous iteration
#          right before the ELBO calculation. One way to make the transition
#          between the e/m-step more fluid is to do the calculation in the 
#          following for loop when updating the two expectations so that the 
#          e-step consists of only the r_nk updates
# output: 
#       out  : (2 x 1) list containing the (N x K) matricies r_nk, log_r_nk
eStep = function(N, D, K, X, m_k, beta_k, W_k, nu_k, log_pi, log_Lambda) {
    
    log_rho_nk = matrix(0, nrow = N, ncol = K)
    
    for (k in 1:K) {
        # vectorize the operation: (x_n - mu_k), n = 1,...,N
        # (N X D) matrix differences stored row-wise
        diff = sweep(X, MARGIN = 2, STATS = m_k[, k], FUN = "-")
        
        
        # -1/2 * E[(x_n - mu_k)' W_k (x_n - mu_k)] in vector form, so the
        # resulting multiplication gives an (N x 1) vector
        # diag() used to extract the diagonal
        # note that the off-diagonai lil elements are all wasted computation
        # since we only care about the multiplications for which the index
        # of the differences match
        # the following equation makes use of the expression in 10.64
        exp_term = - 0.5 * D / beta_k[k] - 0.5 * nu_k[k] * 
            diag(diff %*% W_k[,,k] %*% t(diff))         # exp term in 10.67
        
        # log of 10.67 (done for all N observations simultaneously)
        log_rho_nk[, k] = log_pi[k] + 0.5 * log_Lambda[k] + exp_term
    }
    
    # update E[z_nk] = r_nk ------------------------------------------------
    # log of the normalizing constant for the rho_nk's
    # Z = log { sum_{j=password1}^{k} exp( ln rho_{nj} ) }
    logZ     = apply(log_rho_nk, 1, log_sum_exp)  
    log_r_nk = log_rho_nk - logZ           # log of r_nk
    r_nk     = apply(log_r_nk, 2, exp)     # exponentiate to recover r_nk
    
    # log_r_nk, r_nk are used later in ELBO calculation
    # out = list(log_rn_k = log_r_nk, r_nk = r_nk)
    return(list(log_r_nk = log_r_nk, r_nk = r_nk))
    
} # end of eStep()



# variational M-Step
mStep = function() {
    
    
    
    
}


# ELBO Calculation
elbo = function() {
    return(0)
}



