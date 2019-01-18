

## elbo.R -- covariate-INDEPENDENT case

## calculate the variational lower bound using the current variational params

library(matrixcalc)

elbo = function(theta, prior) {
    
    # initialize the 7 expectations to be computed
    e_ln_p_y = e_ln_p_z = e_ln_p_gamma = e_ln_p_beta_tau = 0
    e_ln_q_z = e_ln_q_beta_tau = e_ln_q_gamma = 0
    
    X = prior$X
    y = prior$y
    N = prior$N
    D = prior$D
    K = prior$K

    # compute 7 expectations ---------------------------------------------------
    
    #### CODE BELOW NEEDS TO BE ADAPTED ----------------------------------------
    
    ln_p_x = ln_p_z = ln_p_pi = ln_p_mu_lambda = 0
    ln_q_z = ln_q_pi = ln_q_mu_lambda = 0
    
    for (k in 1: K) {
        # 10.71 -- 2 * E[ln p(X | Z, mu, Sigma)]
        ln_p_x = ln_p_x + N_k[k] * 
            (log_Lambda[k] - D / beta_k[k] - 
                 nu_k[k] * matrix.trace(S_k[,,k] %*% W_k[,,k]) - 
                 nu_k[k] * t(x_bar_k[, k] - m_k[,k]) %*% W_k[,,k] %*% 
                 (x_bar_k[, k] - m_k[,k]) - D * log(2 * pi))
        
        # 10.74 (partial) -- k-summation terms of the expectation
        ln_p_mu_lambda = ln_p_mu_lambda + D * log(0.5 * beta_0 / pi) + 
            log_Lambda[k] - D * beta_0 / beta_k[k] - 
            beta_0 * nu_k[k] * t(m_k[, k] - m_0) %*% W_k[,,k] %*% 
            (m_k[, k] - m_0) - nu_k[k] * matrix.trace(W_0_inv %*% W_k[,,k])
        
        # 10.77 : E [ ln q(mu, Lambda) ]
        # H(Lambda_k)
        H_k = logB(W = W_k[,,k], nu = nu_k[k]) + 
            0.5 * (nu_k[k] - D - 1) * log_Lambda[k] - 0.5 * nu_k[k] * D
        ln_q_mu_lambda = ln_q_mu_lambda - H_k + 
            0.5 * (log_Lambda[k] + D * log(beta_k[k] / (2 * pi)) - D)
    }
    
    # 10.71  :  E [ ln p(X | Z, mu, Sigma) ] 
    ln_p_x = 0.5 * ln_p_x
    
    # 10.72  :  E [ ln p(Z | pi) ]
    ln_p_z = sum(r_nk %*% log_pi) # sum { (N x K) (K x 1) } --> sum { (N x 1) }
    
    # 10.73  :  E [ ln p(pi) ]
    # same prior concentration parameter for all K clusters
    ln_p_pi = lgamma(K * alpha_0) - K * lgamma(alpha_0) + 
        sum((alpha_0 - 1) * log_pi)                  
    
    # 10.74  :  E [ ln p(mu, Lambda) ]
    ln_p_mu_lambda = 0.5 * (ln_p_mu_lambda + (nu_0 - D - 1) * sum(log_Lambda)) + 
        K * logB(W_0, nu_0) 
    
    # 10.75  :  E [ ln q(Z) ]
    ln_q_z = sum(r_nk * log_r_nk)
    
    # 10.76  :  E [ ln q(pi) ]
    # alpha_k different for each k
    ln_q_pi = lgamma(sum(alpha)) - sum(lgamma(alpha)) + 
        sum((alpha - 1) * log_pi)
    
    # 10.77 : E [ ln q(mu, Lambda) ] -- already calculated in the previous loop
    
    # 10.70
    vlb = ln_p_x + ln_p_z + ln_p_pi + ln_p_mu_lambda - 
        ln_q_z - ln_q_pi - ln_q_mu_lambda                 
    
    #### CODE ABOVE NEEDS TO BE ADAPTED ----------------------------------------
    
    vlb = 0
    
    return(vlb)
    
} # end vb_elbo() function
