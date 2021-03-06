


## computeELBO.R -- compute the variational lower bound for covariate dependent
## mixture model with sparsity assumption

source('misc.R')


elbo_vs = function(prior, theta) {
    
    p1 = p2 = p3 = p4 = p5 = 0
    q1 = q2 = q3 = q4 = 0
    
    X = prior$X
    y = prior$y
    N = prior$N
    D = prior$D
    K = prior$K
    
    r_nk = theta$r_nk
    
    # precompute some commonly computed quantites
    dig_a = digamma(theta$a_k)
    dig_b = digamma(theta$b_k)
    ak_bk = theta$a_k / theta$b_k        # (K x 1)
    X_mu  = X %*% theta$mu_k             # (N x K) : [X*mu1, X*mu2, ..., X*muK]
    
    
    # (K x 1) : diagonal elements -> (mu1'mu1, mu2'mu2, ..., muK'muK)
    mu_mu = diag(t(theta$mu_k) %*% theta$mu_k)
    
    # (D x 1) : diagonal elements of the (D x D) matrix -> (m1'm1, ..., mD'mD)
    # Recall: m_d is (K x D), so m_d'm_d is (D x D)
    md_md = diag(t(theta$m_d) %*% theta$m_d)
    
    # (1) E [ln p (y | - )] =: p1 ----------------------------------------------
    
    p1_tmp0 = numeric(K)
    p1_tmp1 = p1_tmp2 = p1_tmp3 = numeric(D)
    for (k in 1:K) {
        p1_tmp0[k] = sum(r_nk[,k] * (log(2 * pi) - (dig_a - dig_b) + 
                                         ak_bk * prior$y_sq))
    }
    
    for (d in 1:D) {
        p1_tmp1[d] = quadMult(theta$m_d[,d], diag(theta$U_d[,d]))
        p1_tmp2[d] = matrix.trace(diag(theta$U_d[,d]) %*% theta$Q_d_inv[,,d])
        
        j_index = c(1:D)[-d]
        j_sum = numeric(K)
        for (j in j_index) {
            j_sum = j_sum + 
                theta$lambda_d[j] * diag(theta$R_dj[,j]) %*% theta$m_d[,j]
        }
        
        p1_tmp3[d] = t(theta$m_d[,d]) %*% (theta$zeta_d[,d] - 0.5 * j_sum)
    }
    
    p1 = -0.5 * (sum(p1_tmp0) + 
                     sum(theta$lambda_d * (p1_tmp1 + p1_tmp2 + p1_tmp3)))
    
    # (2) E [ln p (Z | - )] =: p2 ----------------------------------------------
    
    p2_tmp = numeric(K)
    for (k in 1:K) {
        # swap summations -> perform summation over n for each k
        p2_tmp[k] = sum(r_nk[,k] * (X_mu[,k] - theta$alpha - theta$phi))
    }
    
    p2 = sum(p2_tmp) # vectorized summation over k
    
    
    # (3) E [ln p (gamma)] =: p3 -----------------------------------------------
    
    p3_tmp = numeric(K)
    for (k in 1:K) {
        p3_tmp[k] = mu_mu[k] + matrix.trace(as.matrix(theta$V_k_inv[,,k]))
    }
    
    p3 = - 0.5 * K * D * log(2 * pi) - 0.5 * sum(p3_tmp)
    
    # (4) E [ln p (tau)] =: p4 -------------------------------------------------
    p4 = K * (prior$a_0 * log(prior$b_0) - lgamma(prior$a_0)) - 
        sum(prior$b_0 * ak_bk + (prior$a_0 - 1) * (dig_a - dig_b))
    
    
    # (5) E [ln p (beta, omega)] =: p5 -----------------------------------------
    p5_tmp = numeric(D)
    for (d in 1:D) {
        p5_tmp[d] = prior$xi_0 * (matrix.trace(as.matrix(theta$Q_d_inv[,,d])) +
                                      md_md[d])
    }
    p5_tmp = theta$lambda_d * log(1 - prior$pi_d) + 
        0.5 * theta$lambda_d * (K * log(prior$xi_0) - K * log(2 * pi) - p5_tmp)
    
    p5 = sum(p5_tmp)
    
    # (6) E [ln q (beta, omega)] =: q1 -----------------------------------------
    q1_tmp = numeric(D)
    for (d in 1:D) {
        q1_tmp[d] = log(det(as.matrix(theta$Q_d[,,d])))
        cat("q1_tmp = ", q1_tmp[d], '\n')
    }
    
    cat("q1_1 = ", theta$lambda_d * log(theta$lambda_d), '\n')
    cat("q1_2 = ", (1 - theta$lambda_d) * log(1 - theta$lambda_d), '\n')
    cat("q1_3 = ", 0.5 * theta$lambda_d * (- K * log(2 * pi) + q1_tmp), '\n')
    
    q1 = sum(theta$lambda_d * log(theta$lambda_d) + 
        (1 - theta$lambda_d) * log(1 - theta$lambda_d) + 
        0.5 * theta$lambda_d * (- K * log(2 * pi) + q1_tmp))
    
    
    
    # (7) E [ln q (gamma)] =: q2 -----------------------------------------------
    q2_tmp = numeric(K)
    for (k in 1:K) {
        q2_tmp[k] = log(det(as.matrix(theta$V_k[,,k])))
    }
    q2 = - 0.5 * K * D * (log(2 * pi) + 1)  + 0.5 * sum(q2_tmp)
    
    # (8) E [ln q (Z)] =: q3 ---------------------------------------------------
    q3 = sum(r_nk * theta$log_r_nk)
    
    # (9) E [ln q (tau)] =: q4 -------------------------------------------------
    q4 = sum(theta$a_k * (theta$b_k - 1) - lgamma(theta$a_k) + 
        (theta$a_k - 1) * (dig_a - dig_b))
    
    
    # print(paste('E [ln p (y | - )] =', p1))
    cat('E [ln p (y | - )]',      '=',   p1,  '\n',  sep = '\t')
    cat('E [ln p (Z | - )]',      '=',   p2,  '\n',  sep = '\t')
    cat('E [ln p (gamma)]',       '=',   p3,  '\n',  sep = '\t')
    cat('E [ln p (tau)]\t\t',     '=\t', p4,  '\n',  sep = '')
    cat('E [ln p (beta, omega)]', '=',   p5,  '\n',  sep = '\t')
    
    cat('E [ln q (beta, omega)]', '=',   q1,  '\n',  sep = '\t')
    cat('E [ln q (gamma)]',       '=',   q2,  '\n',  sep = '\t')
    cat('E [ln q (Z)]\t\t',       '=\t', q3,  '\n',  sep = '')
    cat('E [ln q (tau)]\t\t',     '=\t', q4,  '\n',  sep = '')
    
    elbo = p1 + p2 + p3 + p4 + p5 - q1 - q2 - q3 - q4
    
    return(elbo)
    
    
}





