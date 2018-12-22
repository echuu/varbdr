
library(Matrix)
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
eStep = function(N, D, K, X, theta) {
    
    # initialize local variational parameters
    m_k        = theta$m_k
    beta_k     = theta$beta_k
    W_k        = theta$W_k
    nu_k       = theta$nu_k
    log_pi     = theta$log_pi
    log_Lambda = theta$log_Lambda
    
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
mStep = function(N, K, D, X, theta, alpha_0, beta_0, nu_0, W_0_inv, m0) {
    
    x_bar_k = theta$x_bar_k
    r_nk = theta$r_nk
    m_k = theta$m_k
    W_k = theta$W_k
    W_k_inv = theta$W_k_inv
    S_k = theta$S_k
    log_Lambda = theta$log_Lambda
    log_pi = theta$log_pi
    
    # calculate the quantities: N_k, xbar_k, S_k, used in the updates of
    # the variational distributions
    N_k = colSums(r_nk) + 1e-10  # (K x 1) : N_1, N_2, ... , N_K ----  10.51
    for (k in 1:K) {
        x_bar_k[, k] = (r_nk[ ,k] %*% X) / N_k[k]        # (D x K) -- 10.52
        # x_diff = (x_n - xbar_k) for n = 1,..,N; ---- (N x D)
        x_diff     = sweep(X, MARGIN = 2, STATS = x_bar_k[, k], FUN = "-")
        # vectorized calculation of S_k as shown in 10.53 -- (D x D)
        # store the resulting matrix in the k-th index of the array
        # S_k[, , k] = t(x_diff) %*% diag(r_nk[,k] / N_k[k]) %*% x_diff
        # same calculation as above, uses less memory    # (D x D) -- 10.53
        S_k[,,k] = t(x_diff) %*% (x_diff * r_nk[, k]) / N_k[k]  # 
    }
    
    ## Update Dirichlet parameter
    alpha = alpha_0 + N_k                                # (K x 1) -- 10.58
    ## Expectation of mixing proportions: E(pi_k)            
    pi_k = (alpha_0 + N_k) / (K * alpha_0 + N)           # (K x 1) -- 10.69
    
    ## Update parameters for G-W distribution
    beta_k = beta_0 + N_k                                # (K x 1) -- 10.60
    nu_k   = nu_0 + N_k + 1                              # (K x 1) -- 10.63
    for (k in 1:K) {
        # update mean parameter for beta_k --              (K x 1) -- 10.61
        m_k[, k] = (1 / beta_k[k]) * (beta_0 * m_0 + N_k[k] * x_bar_k[, k])  
        
        # update variance parameter for beta_k --          (D x D) -- 10.62
        W_k_inv = W_0_inv + N_k[k] * S_k[,,k] + 
            ((beta_0 * N_k[k]) / (beta_0 + N_k[k])) * 
            tcrossprod((x_bar_k[, k] - m_0))
        # previous matrix is inverted matrix that we want
        W_k[, , k] = solve(W_k_inv)                    
    }
    
    
    ## Update expectations over \pi and \Lambda
    # E [ log(det(Lambda_k)) ]
    for (k in 1:K) {                                               
        log_Lambda[k] = sum(digamma((nu_k[k] + 1 - 1:D) / 2)) + 
            D * log(2) + log(det(W_k[, , k]))            # (K x 1) -- 10.65
    }
    # E[ log(pi_k) ], k = 1,...,K
    log_pi = digamma(alpha) - digamma(sum(alpha))        # (K x 1) -- 10.66
    
    theta = list()
    
    # store the updated expectations
    theta$log_Lambda = log_Lambda
    theta$log_pi = log_pi
    
    # store defined variables
    theta$N_k = N_k
    theta$S_k = S_k
    theta$x_bar_k = x_bar_k
    
    # store variational parameters
    theta$m_k = m_k
    theta$W_k = W_k
    theta$nu_k = nu_k
    theta$beta_k = beta_k
    theta$alpha = alpha
    theta$pi_k = pi_k
    
    # theta$r_nk = r_nk
    # theta$log_r_nk = theta$log_r_nk
    
    return(theta)
    
} # end of mStep()


# ELBO Calculation
elbo = function() {
    return(0)
}



