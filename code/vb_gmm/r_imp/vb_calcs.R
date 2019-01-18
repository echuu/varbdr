
library(Matrix)
source("misc.R")


# initialize variational parameters + related variables
# intput: 
#    log_lambda  : (K x 1) E [ log(det(Lambda_k)) ]
#    log_pi      : (K x 1) E [ log(pi) ]
#    N_k
#    S_k
#    x_bar_k
#    m_k
#    W_k
#    nu_k
#    beta_k
#    alpha
#    r_nk
#    log_r_nk
initVarParams = function(K, log_Lambda, log_pi, N_k, S_k, x_bar_k, 
                         m_k, W_k, nu_k, beta_k, alpha, r_nk, log_r_nk) {
    
    theta = list()
    # store the updated expectations
    theta$log_Lambda = log_Lambda
    theta$log_pi = log_pi
    
    # store defined variables
    theta$N_k = numeric(K)
    theta$S_k = S_k
    theta$x_bar_k = x_bar_k
    
    # store variational parameters
    theta$m_k = m_k
    theta$W_k = W_k
    theta$nu_k = nu_k
    theta$beta_k = beta_k
    theta$alpha = alpha
    theta$pi_k = numeric(K)
    
    # e-step parameters
    theta$r_nk = r_nk
    theta$log_r_nk = log_r_nk
    
    return(theta)
} # end initVarParams()



# variational E-Step -- update r_nk using the re-estimated variational
#                       parameters from the previous iteration's m-step
# input:
#       N          : (1 x 1) number of observations
#       D          : (1 x 1) dimension of covariates
#       K          : ()
#       X          : (N x D) design matrix
#       theta      : list containing the variational parameters
# output: 
#       theta      : list containing variational parameters with updated
#                    quantity of r_nk, log_r_nk
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
    
    
    # Essential E-Step: update r_nk
    theta$r_nk = r_nk
    theta$log_r_nk = log_r_nk
    
    # log_r_nk, r_nk are used later in ELBO calculation
    return(theta)
    
} # end of eStep()



# variational M-Step
mStep = function(N, K, D, X, theta, alpha_0, beta_0, nu_0, W_0_inv, m_0) {
    
    x_bar_k    = theta$x_bar_k
    r_nk       = theta$r_nk
    m_k        = theta$m_k
    W_k        = theta$W_k
    W_k_inv    = theta$W_k_inv
    S_k        = theta$S_k
    log_Lambda = theta$log_Lambda
    log_pi     = theta$log_pi
    
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
    ## ***** these 2 (of the total 3 expectations) are updated here since
    ## they are used in the calculation of the ELBO
    ## The last expectation, which isn't needed in the ELBO, is computed
    ## as needed in the NEXT iteration in the E-step (I think the purpose of
    ## splitting these calculations is to save the extra computation in case
    ## the elbo converges, save one expectation calculation)
    ## see elbo below to verify
    
    # E [ log(det(Lambda_k)) ]
    for (k in 1:K) {                                               
        log_Lambda[k] = sum(digamma((nu_k[k] + 1 - 1:D) / 2)) + 
            D * log(2) + log(det(W_k[, , k]))            # (K x 1) -- 10.65
    }
    # E[ log(pi_k) ], k = 1,...,K
    log_pi = digamma(alpha) - digamma(sum(alpha))        # (K x 1) -- 10.66
    
    # theta = list()
    
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
    
    return(theta)
    
} # end of mStep()


# elbo(): calculate the variational lower bound
## input: 
    # D          : (1 x 1) dimension of coefficient vector
    # K          : (1 x 1) # of clusters
    # theta      : list containing variational parameters + other model vars
        # r_nk       : (N x K) responsibilities for (n, k)
        # log_r_nk   : (N x K) log of responsibilities
        # N_k        : (K x 1) sum of responsibilities for each cluster k
        # S_k        : K-dim array of  (D x D) matrices as shown in 10.53
        # x_bar_k    : (D x K) weighted mean (indexed by k) stored column-wise
        # W_k        : K-dim array of (D x D) scale matrix for lambda
        # nu_k       : (K x 1) degrees of freedom for each Lambda_k
        # m_k        : (D x K) mean for mu_k stored in the k-th column as D-dim vec
        # beta_k     : (K x 1) variance scale for each beta_k
        # alpha      : (K x 1) concentration parameters
        # log_pi     : (K x 1) E [ log pi ]
        # log_Lambda : (K x 1) E [ log(det(Lambda_k)) ]
    # alpha_0    : (1 x 1) prior concentration parameter (symmetric)
    # m_0        : (D x 1) prior mean for beta
    # beta_0     : (1 x 1) prior variance scale for beta
    # W_0        : (D x D) prior scale matrix for lambda
    # W_0_inv    : (D x D) inverse of prior scale matrix for lambda
    # nu_0       : (1 x 1) prior degrees of freedom for each Lambda_k
## output: 
# elbo    : variational lower bound for current values of variational params
elbo = function(D, K, theta, # r_nk, log_r_nk,
                alpha_0, m_0, beta_0, W_0, W_0_inv, nu_0) {
    
    log_r_nk = theta$log_r_nk
    r_nk     = theta$r_nk
    
    
    # helpful definitions
    N_k     = theta$N_k
    S_k     = theta$S_k
    x_bar_k = theta$x_bar_k
    
    # update variational parameters
    W_k    = theta$W_k 
    nu_k   = theta$nu_k 
    m_k    = theta$m_k
    beta_k = theta$beta_k
    alpha  = theta$alpha
    pi_k   = theta$pi_k
    
    # update expectations
    log_Lambda = theta$log_Lambda
    log_pi     = theta$log_pi
    
    
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
    
    return(vlb)
    
} # end of calculateELBO() function


## various checks on the ELBO
## input:
# VERBOSE      : logical, if TRUE, then progress printed each iter
# i            : current iteration
# max_iter     : max # of iterations
# epsilon_conv : tolerance used to assess convergence
checkELBO = function(VERBOSE, i, max_iter, L, epsilon_conv) {
    
    CONVERGE = FALSE
    
    # display iteration, ELBO, change in ELBO
    if (VERBOSE) {
        cat("It:\t",i,"\tLB:\t",L[i], "\tLB_diff:\t",L[i] - L[i - 1],"\n")
    }
    
    # Check if lower bound decreases
    if (L[i] < L[i - 1]) { 
        message("Warning: Lower bound decreases!\n")
    }
    # Check for convergence
    if (abs(L[i] - L[i - 1]) < epsilon_conv) { 
        CONVERGE = TRUE 
    }
    
    # Check if VB converged in the given maximum iterations
    if (i == max_iter) {
        warning("VB did not converge!\n")
    }
    
    return(CONVERGE)
} # end of checkELBO() function





