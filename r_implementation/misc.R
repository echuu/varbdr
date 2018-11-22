
## misc.R

library(matrixcalc)


## commonly used calculations in Bayesian computation, including some tricks
## to stabilize calculations

# log_sum_exp():
# calculates expressions of the form log(sum(exp(x)))
log_sum_exp = function(x) { 
    offset = max(x)                         # scale by max to prevent overflow
    s = log(sum(exp(x - offset))) + offset
    i = which(!is.finite(s))                # check for any overflow
    if (length(i) > 0) {                    # replace inf values with max
        s[i] = offset 
    }
    
    return(s)
} # end of log_sum_exp()


# Wishart related calculations

## logB() function
# bishop appendix B.79
# log(B(W, nu)), B() is the normalizing constant of the Wishart distribution
logB = function(W, nu) {
    D = ncol(W)
    det_term = - 0.5 * nu * log(det(W))
    power_term = - (0.5 * nu * D * log(2) + 0.25 * D *(D - 1) * log(pi))
    prod_term = - sum(lgamma(0.5 * (nu + 1 - 1:D)))
    return(det_term + power_term + prod_term)
}

## logDirConst()
# log of normalizing constant for dirichlet(alpha)
logDirConst = function(alpha) {
    return(lgamma(sum(alpha)) - sum(lgamma(alpha)))
}


########## ----------    problem-specific calculations    ---------- ###########


# calculateELBO(): calculate the variational lower bound
## input: 
    # D          : (1 x 1) dimension of coefficient vector
    # K          : (1 x 1) # of clusters
    # r_nk       : (N x K) responsibilities for (n, k)
    # log_r_nk   : (N x K) log of responsibilities
    # N_k        : (K x 1) sum of responsibilities for each cluster k
    # S_k        : K-dim array of  (D x D) matrices as shown in 10.53
    # x_bar_k    : (D x K) weighted mean (indexed by k) stored column-wise
    # alpha_0    : (1 x 1) prior concentration parameter (symmetric)
    # m_0        : (D x 1) prior mean for beta
    # beta_0     : (1 x 1) prior variance scale for beta
    # W_0        : (D x D) prior scale matrix for lambda
    # W_0_inv    : (D x D) inverse of prior scale matrix for lambda
    # nu_0       : (1 x 1) prior degrees of freedom for each Lambda_k
    # W_k        : K-dim array of (D x D) scale matrix for lambda
    # nu_k       : (K x 1) degrees of freedom for each Lambda_k
    # m_k        : (D x K) mean for mu_k stored in the k-th column as D-dim vec
    # beta_k     : (K x 1) variance scale for each beta_k
    # alpha      : (K x 1) concentration parameters
    # log_pi     : (K x 1) E [ log pi ]
    # log_Lambda : (K x 1) E [ log(det(Lambda_k)) ]
## output: 
    # elbo    : variational lower bound for current values of variational params
calculateELBO = function(D, K, r_nk, log_r_nk, N_k, S_k, x_bar_k, alpha_0, 
                         m_0, beta_0, W_0, W_0_inv, nu_0, W_k, nu_k, 
                         m_k, beta_k, alpha, log_pi, log_Lambda) {
    
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
    elbo = ln_p_x + ln_p_z + ln_p_pi + ln_p_mu_lambda - 
        ln_q_z - ln_q_pi - ln_q_mu_lambda                       
    
    return(elbo)
    
} # end of calculateELBO() function




