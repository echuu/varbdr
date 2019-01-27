
## eStep.R
## perform the variational e-step -- covariate-INDEPENDENT case

library(matrixcalc)
source("misc.R")

# input: 
#         theta : list of variational parameters
# output:
#         theta : list of variational parameters with 
#                 r_nk, log_r_nk, updated
eStep = function(theta, prior) {
    
    # print("e-step")
    
    X = prior$X
    y = prior$y
    N = prior$N
    K = prior$K
    D = prior$D
    
    # initialize local copies of variational parameters to be updated:
    r_nk       = matrix(0, N, K)       # (N x K) : normalized responsibilities
    log_r_nk   = matrix(0, N, K)       # (N x K) : log(r_nk)
    log_rho_nk = matrix(0, N, K)       # (N x K) : log(rho_nk)
    
    # these 3 quantities can be potentially saved in theta so they need not
    # be recalculated every time we need them (fairly common)
    # if we want to make this more efficient, we should: 
        # (1) calculate once in initiVarParams() so that the first iteration 
        #     of CAVI has access to these values
        # (2) calculate and store in theta during mStep() so that the NEXT iter
        #     of cAVI has access to these values
    # do this after the algorithm is up and running
    psi_a = digamma(theta$a_k)         # (K x 1)
    psi_b = digamma(theta$b_k)         # (K x 1)
    ak_bk = theta$a_k / theta$b_k      # (K x 1) : a_k / b_k
    
    # E[ln pi_k] --- vectorized calculation for pi_1, .., pi_K
    e_ln_pi = digamma(theta$alpha_k) - digamma(sum(theta$alpha_k)) # (K x 1)
    
    # X_mu  = X %*% theta$m_k            # (N x K) : (N x D) * (D x K)
    
    # update log_rho_nk (row-wise updates)
    for (n in 1:N) {
        x_Vinv_x = numeric(K)
        for (k in 1:K) {
            # print(length(X[n,]))
            # print(dim(theta$V_k_inv[,,k]))
            x_Vinv_x[k] = t(X[n,]) %*% theta$V_k_inv[,,k] %*% X[n,] # scalar
        } # end of inner for()
        
        # populate the n-th row with a k-dim vector
        # first line in calculation below is the only difference from the
        # covariate-dependent case
        log_rho_nk[n,] = e_ln_pi -
            (log(2 * pi) - psi_a + psi_b + x_Vinv_x + 
                    ak_bk * (y[n] - t(X[n,] %*% theta$m_k))^2) / 2
        
    } # end of outer for()
    
    for (n in 1:N) {
        for (k in 1:K) {
            
            log_rho_nk[n,k] = -0.5 * log(2 * pi) + 0.5 * (psi_a[k] - psi_b[k]) -
                0.5 * (theta$a_k[k] / theta$b_k[k] * 
                           (y[n] - t(X[n,]) %*% theta$m[,k]) + 
                           t(X[n,]) %*% theta$V_k_inv[,,k] %*% X[n,]) +
                digamma(theta$alpha_k[k]) - digamma(sum(theta$alpha_k))
            
        }
    }
    
    # rho_nk = exp(log_rho_nk)
    
    # compute r_nk = divide each element of rho_nk by the sum of the
    #                corresponding row; each row consists of 'responsibilities'
    #                of each of the clusters for the observation in that row
    # r_nk = sweep(rho_nk, MARGIN = 1, STAT = rowSums(rho_nk), FUN = '/')
    
    
    # compute log_r_nk
    # log_r_nk = log(r_nk)
    
    
    logZ     = apply(log_rho_nk, 1, log_sum_exp)  # log of normalizing constant
    log_r_nk = log_rho_nk - logZ                  # log of r_nk
    r_nk     = apply(log_r_nk, 2, exp)            # exponentiate to recover r_nk
    
    # cat("difference in r_nk:", sum(theta$r_nk - r_nk), "\n")
    
    theta$log_r_nk    = log_r_nk                # (N x K)
    theta$r_nk        = r_nk                    # (N x K)
    theta$N_k         = colSums(r_nk) + 1e-10   # (K x 1)
    
    return(theta)
    
} # end eStep()


# end of eStep.R

