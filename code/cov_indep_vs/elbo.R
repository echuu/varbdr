
# elbo.R --- implementation of the calculation of the variational lower bound



# computeELBO() : compute the variational lower bound using the current settings
#                 of the variational family, compare value of lower bound w/
#                 previous iteration's lower bound and determine convergence
#                 status
computeELBO = function(prior, theta) {
    
    # ELBO is a function of the following 9 expectations: 
        # (1) E [ln p (y | - )]      =: p1
        # (2) E [ln p (Z | - )]      =: p2
        # (3) E [ln p (pi)]          =: p3
        # (4) E [ln p (tau)]         =: p4
        # (5) E [ln p (beta, omega)] =: p5
        # (6) E [ln q (beta, omega)] =: q1
        # (7) E [ln q (pi)]          =: q2
        # (8) E [ln q (Z)]           =: q3
        # (9) E [ln q (tau)]         =: q4
    
    p1 = p2 = p3 = p4 = p5 = 0
    q1 = q2 = q3 = q4 = 0
    
    X = prior$X
    y = prior$y
    N = prior$N
    D = prior$D
    K = prior$K  
    
    r_nk = theta$r_nk
    l_d  = theta$lambda_d
    X_mk  = X %*% theta$m_k        # (N x K)
    
    # precompute some commonly computed quantites
    dig_a = digamma(theta$a_k)
    dig_b = digamma(theta$b_k)
    ak_bk = theta$a_k / theta$b_k        # (K x 1)
    
    # E_q [ln pi_k]
    # e_ln_pi = sum(digamma(theta$alpha_k) - digamma(sum(theta$alpha_k)))
    
    # (K x 1) : expectation of each of the pi_k's
    e_ln_pi_k = digamma(theta$alpha_k) - digamma(sum(theta$alpha_k))
    
    
    # (D x 1) : diagonal elements of the (D x D) matrix -> (m1'm1, ..., mD'mD)
    # Recall: m_d is (K x D), so m_d'm_d is (D x D)
    md_md = diag(t(theta$m_d) %*% theta$m_d)
    
    # (1) E [ln p (y | - )] =: p1 ----------------------------------------------
    
    # p1_tmp0 = numeric(K)
    # p1_tmp1 = p1_tmp2 = p1_tmp3 = numeric(D)
    # for (k in 1:K) {
    #     p1_tmp0[k] = sum(r_nk[,k] * (log(2 * pi) - (dig_a - dig_b) + 
    #                                      ak_bk * prior$y_sq))
    # }
    
    # for (d in 1:D) {
    #    p1_tmp1[d] = quadMult(theta$m_d[,d], diag(theta$U_d[,d]))
    #    p1_tmp2[d] = matrix.trace(diag(theta$U_d[,d]) %*% theta$Q_d_inv[,,d])
        
    #    j_index = c(1:D)[-d]
    #    j_sum = numeric(K)
    #    for (j in j_index) {
    #        j_sum = j_sum + 
    #            theta$lambda_d[j] * diag(theta$R_dj[,j]) %*% theta$m_d[,j]
    #   }
        
    #    p1_tmp3[d] = t(theta$m_d[,d]) %*% (theta$zeta_d[,d] - 0.5 * j_sum)
    # }
    
    # p1 = -0.5 * (sum(p1_tmp0) + 
    #                  sum(theta$lambda_d * (p1_tmp1 + p1_tmp2 + p1_tmp3)))
    
    e1 = numeric(K)
    for (k in 1:K) {
        
        xSigmax = numeric(N)
        for (n in 1:N) {
            xSigmax[n] = t(X[n,]) %*% theta$Sigma_k[,,k] %*% X[n,]
        }
        
        # all quantities in the sum should be 1-dim or N-dim
        e1[k] = sum(r_nk[,k] * (log(2 * pi) - (dig_a[k] - dig_b[k]) + 
                                    ak_bk[k] * ((y - X_mk[,k])^2 + xSigmax)))
    }
    e_ln_p_y = -0.5 * sum(e1) # outer sum over k
    
    
    # (2) E [ln p (Z | - )] = --------------------------------------------------
    # p2_tmp = numeric(N) 
    # for (n in 1:N) {
    #     p2_tmp[n] = sum(r_nk[n,] * e_ln_pi) # summation over k
    # }
    # p2 = sum(p2_tmp) # summation over n
    
    e2 = numeric(K)
    for (k in 1:K) {
        e2[k] = sum(r_nk[,k] * (e_ln_pi_k[k])) # summation over n
    }
    
    e_ln_p_z = sum(e2) # outer sum over k
    
    
    # (3) E [ln p (pi)] = ------------------------------------------------------
    log_C_alpha_0 = log_dir_const(prior$alpha_0, K)
    e_ln_p_pi = log_C_alpha_0 + sum((prior$alpha_0 - 1) * e_ln_pi)
    
    
    # (4) E [ln p (beta, omega)] = ---------------------------------------------
    
    # p5_tmp = numeric(D)
    # for (d in 1:D) {
    #    p5_tmp[d] = prior$xi_0 * 
    #        (matrix.trace(as.matrix(theta$Q_d_inv[,,d])) + md_md[d])
    # }
    # p5_tmp = theta$lambda_d * log(1 - prior$pi_d) + 
    #     0.5 * theta$lambda_d * (K * log(prior$xi_0) - K * log(2 * pi) - p5_tmp)
    
    # p5 = sum(p5_tmp)
    
    e4_0 = sum(l_d * log(prior$pi_d) + (1 - l_d) * log(1 - prior$pi_d)) # R^1
    qd_trace = numeric(D)
    for (d in 1:D) {
        qd_trace[d] = matrix.trace(as.matrix(theta$Q_d_inv[,,d]))
    }
    e_ln_p_beta = e4_0 + sum(prior$xi_0 * (qd_trace + md_md[d]))
    
    
    # (5) E [ln p (tau)] = -----------------------------------------------------
    e_ln_p_tau = K * (prior$a_0 * log(prior$b_0) - lgamma(prior$a_0)) - 
        sum(prior$b_0 * ak_bk + (prior$a_0 - 1) * (dig_a - dig_b))
    
    
    # (6) E [ln q (beta, omega)] =: q1 -----------------------------------------
    
    # q1_tmp = numeric(D)
    # for (d in 1:D) {
    #    q1_tmp[d] = log(det(as.matrix(theta$Q_d[,,d])))
    #    cat("q1_tmp = ", q1_tmp[d], '\n')
    # }
    
    # cat("q1_1 = ", theta$lambda_d * log(theta$lambda_d), '\n')
    # cat("q1_2 = ", (1 - theta$lambda_d) * log(1 - theta$lambda_d), '\n')
    # cat("q1_3 = ", 0.5 * theta$lambda_d * (- K * log(2 * pi) + q1_tmp), '\n')
    
    # missing (- K) term in the last summation
    # q1 = sum(theta$lambda_d * log(theta$lambda_d) + 
    #              (1 - theta$lambda_d) * log(1 - theta$lambda_d) + 
    #              0.5 * theta$lambda_d * (- K * log(2 * pi) + q1_tmp))
    
    log_det_qd = numeric(D)
    for (d in 1:D) {
        log_det_qd[d] = log(det(as.matrix(theta$Q_d[,,d])))
    }
    
    # check values of lambda_d for each d
    if (theta$VERBOSE) {
        for (d in 1:D) {
            cat("p(lambda_", d, ') \t\t= \t', l_d[d], '\n', sep = '')
        }
    }
    
    e_ln_q_beta = sum(l_d * log(l_d) + (1 - l_d) * log(1 - l_d) + 
                          0.5 * l_d * (-K * log(2 * pi) + log_det_qd - K))
    
    
    # (7) E [ln q (Z)] = -------------------------------------------------
    e_ln_q_z = sum(r_nk * theta$log_r_nk) # sum over (N x K) matrix
    
    
    # (8) E [ln q (tau)] = -----------------------------------------------------
    e_ln_q_tau = sum(theta$a_k * (log(theta$b_k) - 1) - lgamma(theta$a_k) + 
                         (theta$a_k - 1) * (dig_a - dig_b))
    
    
    # TODO: (9) E [ln q (pi)] = --------------------------------------------------
    log_C_alpha_k = log_dir_const(theta$alpha_k, K)
    
    # previous update is wrong - e_ln_pi already sums over the alpha_k's
    # q2 = log_C_alpha_k + sum((theta$alpha_k - 1) * e_ln_pi) # k-sum
    
    e_ln_q_pi = log_C_alpha_k + sum((theta$alpha_k - 1) * e_ln_pi_k)
    
    
    
    if (theta$VERBOSE) {
        cat('E [ln p (y | - )]',      '=',   e_ln_p_y,     '\n',  sep = '\t')
        cat('E [ln p (Z | - )]',      '=',   e_ln_p_z,     '\n',  sep = '\t')
        cat('E [ln p (pi)]\t',        '=',   e_ln_p_pi,    '\n',  sep = '\t')
        cat('E [ln p (beta, omega)]', '=',   e_ln_p_beta,  '\n',  sep = '\t')
        cat('E [ln p (tau)]\t\t',     '=\t', e_ln_p_tau,   '\n',  sep = '')
        
        cat('E [ln q (beta, omega)]', '=',   e_ln_q_beta,  '\n',  sep = '\t')
        cat('E [ln q (Z)]\t',         '=',   e_ln_q_z,     '\n',  sep = '\t')
        cat('E [ln q (tau)]\t\t',     '=\t', e_ln_q_tau,   '\n',  sep = '')
        cat('E [ln q (pi)]\t\t',      '=\t', e_ln_q_pi,    '\n',  sep = '')
    }
    
    
    elbo = e_ln_p_y +  e_ln_p_z + e_ln_p_pi + e_ln_p_beta + e_ln_p_tau -
        e_ln_q_beta - e_ln_q_z - e_ln_q_tau - e_ln_q_pi
    
    theta$L[theta$curr] = elbo
    
    # determine convergence status
    if (checkELBO(prior, theta)) {
        theta$converge = TRUE
    }
    
    return(theta)
    
} # end computeELBO() function



# end of computeELBO.R 

