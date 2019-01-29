
## mStep.R
## perform the variational m-step -- covariate-INDEPENDENT case

library(matrixcalc)

# input: 
#          theta : list of variational parameters
# output:
#          theta : list of variational parameters with 
#                  variational parameters updated

mStep = function(theta, prior) {
    
    
    PRINT_PARAMS = FALSE
    
    # print("m-step")
    
    X = prior$X
    y = prior$y
    N = prior$N
    K = prior$K
    D = prior$D
    
    I_D        = diag(1, D)                      # (D X D) : identity matrix
    # X_m        = X %*% theta$m_k                 # (N x K) : (N x D) * (D x K)
    # M          = theta$mu_k %*% t(theta$lambda)  # (D x N) : (D x K) * (K x N)
    Lambda0_m0 = prior$Lambda_0 %*% prior$m_0    # (D x 1) : Lambda_0 * m_0    
    
    # update q(pi) = Dir( pi | alpha_k )
    #    parameter to update: alpha_k
    theta$alpha_k = prior$alpha_0 + theta$N_k
    
    # pi_k = posterior mean of of pi_k
    # theta$pi_k = (prior$alpha_0 + theta$N_k) / (K * prior$alpha_0 + N) # (K x 1)
    
    # print(round(theta$pi_k, 4))
    # cat("Sum of mixing parameters =", sum(theta$pi_k), '\n')
    
    # update q(beta_k | tau_k) = N ( beta_k | m_k, (tau_k V_k)^{-1} ) 
    #     parameters to update: V_k, V_k^{-1}, zeta_k, m_k
    for (k in 1:K) {
        
        r_nk_xx = 0
        # r_nk_yx = 0
        for (n in 1:N) {
            r_nk_xx = r_nk_xx + theta$r_nk[n,k] * crossprod(t(X[n,]))
            # r_nk_yx = r_nk_yx + theta$r_nk[n,k] * y[n] * X[n,]
        }

        # theta$V_k[,,k] = prior$Lambda_0 + t(X) %*% diag(theta$r_nk[,k]) %*% X
        
        theta$V_k[,,k]     = prior$Lambda_0 + r_nk_xx
        theta$V_k_inv[,,k] = solve(theta$V_k[,,k])
        theta$zeta_k[,k]   = Lambda0_m0 + t(X) %*% (theta$r_nk[,k] * y)
        theta$m_k[,k]      = theta$V_k_inv[,,k] %*% theta$zeta_k[,k]
    }
    
    
    # update q(tau): a_k, b_k
    theta$a_k = prior$a_0 + 0.5 * theta$N_k    # (K x 1)
    for (k in 1:K) {
        theta$b_k[k] = - t(theta$zeta_k[,k]) %*% theta$V_k_inv[,,k] %*% 
            theta$zeta_k[,k] + sum(theta$r_nk[,k] * y^2)
    }
    
    # this line throws a warning message: vector + array of length 1
    # operation is deprecated
    # theta$b_k = prior$b_0 + 
    #    0.5 * (theta$b_k + t(prior$m_0) %*% prior$Lambda_0 %*% prior$m_0)
    
    theta$b_k = prior$b_0 + 
        0.5 * (theta$b_k + c(t(prior$m_0) %*% prior$Lambda_0 %*% prior$m_0))
    
    
    # update the random variables using posterior means ------------------------
    
    # mixture weights
    theta$pi_k = (prior$alpha_0 + theta$N_k) / (K * prior$alpha_0 + N) # (K x 1)
    
    # coefficient vector
    theta$beta_k = theta$m_k
    
    # precision componenets (posterior mean of tau_k)
    theta$tau_k = theta$a_k / theta$b_k
    
    
    # print values of b_k
    if (PRINT_PARAMS) {
        
        # m_k, V_k are difficult to print in a coherent way, if necessary
        # we can look at alternatives later to view these per iteration
        
        alpha_k_str = vecToStr(theta$alpha_k)
        a_k_str     = vecToStr(theta$a_k)
        b_k_str     = vecToStr(theta$b_k)
        
        cat("alpha_k =\t", alpha_k_str, "\n", sep = '')
        cat("a_k =\t\t",     a_k_str,     "\n", sep = '')
        cat("b_k =\t\t",     b_k_str,     "\n", sep = '')
    }
    
    # update the current iteration
    theta$curr = theta$curr + 1
    
    return(theta)
    
} # end mStep()


# end of mStep.R




