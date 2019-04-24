

## elbo.R -- covariate-INDEPENDENT case

## calculate the variational lower bound using the current variational params
## note: most of the updates here should be very similar, if not identical,
##       to the updates in the covariate-DEPENDENT case, since the cov-dependent
##       case introduces dependency via another coefficient vector (gamma_k's)
##       so those variational updates should have an effect that is more local,
##       in the sense that the updates shouldn't strongly influence other
##       variational updates.

source(paste(COV_INDEP, MISC_FUNCS, sep = '/'))

elbo = function(theta, prior) {
    
    
    PRINT_ELBO_VALUES = FALSE
    
    # initialize the 7 expectations to be computed
    e_ln_p_y = e_ln_p_z = e_ln_p_gamma = e_ln_p_beta_tau = 0
    e_ln_q_z = e_ln_q_beta_tau = e_ln_q_gamma = 0
    
    X = prior$X
    y = prior$y
    N = prior$N
    D = prior$D
    K = prior$K

    # pre-compute commonly computed terms
    psi_a = digamma(theta$a_k)    # (K x 1)
    psi_b = digamma(theta$b_k)    # (K x 1)
    ak_bk = theta$a_k / theta$b_k # (K x 1) : a_k / b_k, elementwise division
    
    
    # expectation of log dirichlet 
    #    E_q [ln pi_k] = digamma(alpha_k) - digamma(sum(alpha)), for all k
    e_ln_pi = digamma(theta$alpha_k) - digamma(sum(theta$alpha_k)) # (K x 1)
    
    
    # compute 7 expectations ---------------------------------------------------
    
    # initialize the 7 expectations to be computed
    e_ln_p_y = e_ln_p_z = e_ln_p_beta_tau = e_ln_p_pi = 0
    e_ln_q_z = e_ln_q_pi = e_ln_q_beta_tau = 0
    
    
    # E[ln p(y | X, beta, tau, Z)] --- same as covariate-dependent case --------
    e1 = numeric(N)
    for (n in 1:N) {
        x_Vinv_x = numeric(K)
        for (k in 1:K) {
            x_Vinv_x[k] = t(X[n,]) %*% theta$V_k_inv[,,k] %*% X[n,] # scalar
        } # end of inner for()
        
        # perform the inner k-summation
        e1[n] = sum(theta$r_nk[n,] * 
                        (log(2 * pi) - psi_a + psi_b +
                             ak_bk * (y[n] - t(X[n,] %*% theta$m_k))^2 +
                             x_Vinv_x))
    }
    # outer n-summation
    e_ln_p_y = -0.5 * sum(e1)
    
    
    # E[ln p(z | pi)] --- different from covariate-dependent case --------------
    e2 = matrix(0, N, K) # (N x K)
    for (n in 1:N) {
        e2[n,] = theta$r_nk[n,] * e_ln_pi # (K x 1) element-wise multiplication
    }
    e_ln_p_z = sum(e2) # sum over entire matrix
    
    
    # E[ln p(beta, tau)] --- same as covariate-dependent case ------------------
    e3 = sum((prior$a_0 + 0.5 * D - 1) * (psi_a - psi_b)) -
        K * (0.5 * D * log(2 * pi) - log(det(prior$Lambda_0)) - 
                 prior$a_0 * log(prior$b_0) + lgamma(prior$a_0))
    
    diff = sweep(theta$m_k, MARGIN = 1, STATS = prior$m_0, FUN = '-')  # (D x K)
    e3_diff_k  = numeric(K)  # store (m_k - m_0)' Lambda_0 (m_k - m_0)
    e3_trace_k = numeric(K)  # store trace(Lambda_0 V_K^{-1})
    for (k in 1:K) {
        e3_diff_k[k]  = quadMult(diff[,k], prior$Lambda_0)
        e3_trace_k[k] = matrix.trace(prior$Lambda_0 %*% theta$V_k_inv[,,k])
    }
    e3_diff_k = ak_bk * (e3_diff_k + prior$b_0) # element-wise multiplication
    
    e_ln_p_beta_tau = e3 - 0.5 * sum(e3_diff_k + e3_trace_k)
    
            
    # E[ln p(pi)] --- does not exist in covariate-dependent case ---------------
    #                 note distinction between E[ln p(pi)] and E[ln pi]
    # compute log of dirichlet constant for pi ~ Dir(alpha_0) (prior)
    log_C_alpha_0 = log_dir_const(prior$alpha_0, 1) # log_dir_const() in misc.R
    e_ln_p_pi = log_C_alpha_0 + sum((prior$alpha_0 - 1) * e_ln_pi)
    

    # E[ln q(z)] --- same as covariate-dependent case --------------------------
    e_ln_q_z = sum(theta$r_nk * theta$log_r_nk) # sum over (N x K) matrix
    
    
    # E[ln q(pi)] --- does not exist in covariate-dependent case ---------------
    # compute log of dirichlet constant for pi ~ Dir(alpha_k) (var. dist.)
    log_C_alpha_k = log_dir_const(theta$alpha_k, K)
    e_ln_q_pi = log_C_alpha_k + sum((theta$alpha_k - 1) * e_ln_pi) # k-sum
    
    
    # E[ln q(beta, tau)] --- same as covariate-dependent case ------------------
    e7_k = (0.5 * D + theta$a_k - 1) * (psi_a - psi_b) + 
        theta$a_k * log(theta$b_k) - theta$a_k - lgamma(theta$a_k)  # (K x 1)
    for (k in 1:K) {
        e7_k[k] = e7_k[k] + log(det(theta$V_k[,,k]))
    }
    e_ln_q_beta_tau = sum(e7_k) - 0.5 * K * D * (log(2 * pi) + 1)

    
    # compute the ELBO, sum of the previous 7 expectations ---------------------
    
    # L = E[ln p(y | X, beta, tau, Z)] + E[ln p(z | pi)] + E[ln p(beta, tau)] + 
    #     E[ln p(pi)] - E[ln q(z)] - E[ln q(pi)] - E[ln q(beta, tau)]
    
    if (PRINT_ELBO_VALUES) {
        cat("len. of E[ln p (y | - )]: ",     (e_ln_p_y), "\n")
        cat("len. of E[ln p (Z | - )]: ",     (e_ln_p_z), "\n")
        cat("len. of E[ln p (beta, tau)]: ",  (e_ln_p_beta_tau), "\n")
        cat("len. of E[ln p (pi)]: ",         (e_ln_p_pi), "\n")
        cat("len. of E[ln q (Z)]: ",          (e_ln_q_z), "\n")
        cat("len. of E[ln q (pi)]: ",         (e_ln_q_pi), "\n")
        cat("len. of E[ln q (beta, tau)]: ",  (e_ln_q_beta_tau), "\n")
    }
    
    vlb = e_ln_p_y + e_ln_p_z + e_ln_p_beta_tau + e_ln_p_pi -
        e_ln_q_z - e_ln_q_pi - e_ln_q_beta_tau
    # cat("ELBO:", vlb, "\n", sep = " ")
    
    return(vlb)
    
} # end elbo() function
