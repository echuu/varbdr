


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
    ak_bk = theta$a_k / theta$b_k      # (K x 1)
    X_mu = X %*% theta$mu_k            # (N x K) : [X*mu_1, X*mu_2, ..., X*mu_K]
    
    # (1) E [ln p (y | - )]
    
    p1_tmp0 = numeric(K)
    p1_tmp1 = p1_tmp2 = p1_tmp3 = numeric(D)
    for (k in 1:K) {
        p1_tmp0[k] = sum(r_nk[,k] * (log(2 * pi) - (dig_a - dig_b) + 
                                         ak_bk * prior$y_sq))
    }
    
    for (d in 1:D) {
        p1_tmp1[d] = quadMult(theta$m_d[,d], theta$U_d[,,d])
        p1_tmp2[d] = matrix.trace(theta$U_d[,,d] %*% theta$Q_d_inv[,,d])
        
        j_index = c(1:D)[-d]
        j_sum = numeric(K)
        for (j in j_index) {
            j_sum = j_sum + lambda_d[j] * theta$R_dj[,,j] %*% theta$m_d[,j]
        }
        
        p1_tmp3[d] = t(theta$m_d[,d]) %*% (theta$zeta_d[,d] - 0.5 * j_sum)
    }
    
    p1 = -0.5 * (sum(p1_tmp0) + 
                     sum(theta$lambda_d * (p1_tmp1 + p1_tmp2 + p1_tmp3)))
    
    # (2) E [ln p (Z | - )]
    
    # (3) E [ln p (gamma)]
    
    # (4) E [ln p (tau)]
    
    # (5) E [ln p (beta, omega)]
    
    # (6) E [ln q (beta, omega)]
    
    # (7) E [ln q (gamma)]
    
    # (8) E [ln q (Z)]
    
    # (9) E [ln q (tau)]
    
    
    elbo = p1 + p2 + p3 + p4 + p5 - q1 - q2 - q3 - q4
    
    return(elbo)
    
    
}





