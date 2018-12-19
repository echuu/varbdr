library(matrixcalc)
library(ggplot2)
library(data.table)
library(purrr)        # overlapping transpose function from data.table
library(mvtnorm)
library(Matrix)


# setwd("~/varbdr/r_implementation")
setwd("C:/Users/chuu/varbdr/cov_independent_code/r_implementation")

source("densityCalculations.R")
source("misc.R")



############################# begin algorithm ##################################

## vb_gmm() function  ----------------------------------------------------------
## input:
    # X       : (N x D) -- observations stored row-wise
    # K       : # of clusters
    # alpha_0 : prior concentration parameter for mixing weights (symmetric)
    # m_0     : prior mean for beta
    # beta_0  : prior variance scale for beta
    # W_0     : (D x D) prior scale matrix for lambda
    # nu_0    : prior degrees of freedom for lambda
    # VERBOSE : if TRUE, then convergence output is shown
## output: 
    # X       : (N x D)
    # K       : # of clusters
    # N       : # of observations
    # D       : dimension of covariates
    # pi_k    : updated mixing weights (K x 1)
    # alpha   : updated concetration parameters
    # r_nk    : updated responsibilities
    # m       : updated mean for beta
    # beta    : updated variance scale for beta
    # W       : K-dim array of (D x D) updated scale matrix for each Lambda_k
    # nu      : updated degrees of freedom for each Lambda_k
    # L       : variational lower bound at each iteration (ELBO)
    # dt_all  : graphic parameter
vb_gmm = function(X, K = 3, alpha_0 = 1 / K, m_0 = c(colMeans(X)), beta_0 = 1, 
                  nu_0 = NCOL(X) + 50, W_0 = diag(100, NCOL(X)), max_iter = 500, 
                  epsilon_conv = 1e-4, is_animation = FALSE, VERBOSE = FALSE) {

    ## animation details -------------------------------------------------------
    if (is_animation) {
        # Variables needed for plotting
        xgrid = seq(from = min(X[,1]) - 2, to = max(X[,1]) + 2, length.out = 80)
        ygrid = seq(from = min(X[,2]) - 8, to = max(X[,2]) + 2, length.out = 80)
        dt = data.table(expand.grid(x = xgrid, y = ygrid))
    }
    dt_all = data.table(x = numeric(), y = numeric(), z = numeric(), 
                        iter = numeric())
    ## animation details -------------------------------------------------------
    
    
    X = as.matrix(X)         # N X D design matrix of covariates
    D = NCOL(X)              # Number of features
    N = NROW(X)              # Number of observations
    L = rep(-Inf, max_iter)  # Store the variational lower bounds
    

    W_0_inv = solve(W_0)     # Compute W_0^{-1} (used in updates for W_k)
    
    # initialize storage matrices/vectors for var. parameters: z_nk, pi_k
    r_nk      = log_r_nk = log_rho_nk = matrix(0, nrow = N, ncol = K) # N x K
    x_bar_k   = matrix(0, nrow = D, ncol = K)         # (D x K) Bishop 10.52
    # S_k is a K-dim array that stores (D x D) (covariance) matrices
    S_k       = W_k = array(0, c(D, D, K))           # (D x D) Bishop 10.53
    log_pi    = log_Lambda = rep(0, K) 
    
    # initialize variational parameters: for mu_k, lambda_k, alpha_k, pi_k
    m_k       = t(kmeans(X, K, nstart = 25)$centers)  # Mean of Gaussian
    beta_k    = rep(beta_0, K)                        # Scale of precision mat
    nu_k      = rep(nu_0, K)                          # Degrees of freedom
    alpha     = rep(alpha_0, K)                       # Dirichlet parameter
    
    # E[log(pi)] -- (K x 1) ---> needed for r_nk update
    log_pi    = digamma(alpha) - digamma(sum(alpha))  

    for (k in 1:K) {
        W_k[,,k] =  W_0  # Scale matrix for Wishart
        
        # E [ log(det(Lambda_k)) ] ---> needed for r_nk update
        log_Lambda[k] = sum(digamma((nu_k[k] + 1 - c(1:D)) / 2)) + 
            D * log(2) + log(det(W_k[ , , k]))
    }
    
    #### -----------------------------------------------------------------------
    
    ## animation details -------------------------------------------------------
    if (is_animation) { # Create animation for initial assignments
        my_z = mixture_pdf_t(model = list(m = m_k, W = W_k, beta = beta_k, 
                                     nu = nu_k, alpha = rep(1/K, K)), data = dt)
        dt_all = rbind(dt_all, dt[, z := my_z] %>% .[, iter := 0])
    }
    ## animation details -------------------------------------------------------
    
    
    # Iterate to find optimal parameters
    for (i in 2:max_iter) {
        
        # Start Variational E-Step ---------------------------------------------
        
        # update responsibilities: involves E[log(det(Lambda_k))], 
        # E[log(det(pi_k))], E[(x_n - mu_k)' W_k (x_n - mu_k)]
        
            # question: do we actually "update" all those expectations here?
            # it looks like the only quantity being updated here is r_nk
            # which involes the 'pre'-set values for the two other expectations
            # answer: one of the remaining expectations is updated in the 
            # following loop, and the rest are (again) updated at end of cavi
        
        for (k in 1:K) {
            # vectorize the operation: (x_n - mu_k), n = 1,...,N
            # (N X D) matrix differences stored row-wise
            diff = sweep(X, MARGIN = 2, STATS = m_k[, k], FUN = "-")
            
            
            # -1/2 * E[(x_n - mu_k)' W_k (x_n - mu_k)] in vector form, so the
            # resulting multiplication gives an (N x 1) vector
              # diag() used to extract the diagonal
              # note that the off-diagonal elements are all wasted computation
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
        
        
        # Finish Variational E-Step --------------------------------------------
        
        # Start Variational M-Step ---------------------------------------------
        
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
        
        
        ### testing this chunk -------------------------------------------------
        ## Update expectation r_nk
        for (k in 1:K) {
            # vectorize the operation: (x_n - mu_k), n = 1,...,N
            # (N X D) matrix differences stored row-wise
            diff = sweep(X, MARGIN = 2, STATS = m_k[, k], FUN = "-")
            
            
            # -1/2 * E[(x_n - mu_k)' W_k (x_n - mu_k)] in vector form, so the
            # resulting multiplication gives an (N x 1) vector
            # diag() used to extract the diagonal
            # note that the off-diagonal elements are all wasted computation
            # since we only care about the multiplications for which the index
            # of the differences match
            # the following equation makes use of the expression in 10.64
            exp_term = - 0.5 * D / beta_k[k] - 0.5 * nu_k[k] * 
                diag(diff %*% W_k[,,k] %*% t(diff))         # exp term in 10.67
            what el
            # log of 10.67 (done for all N observations simultaneously)
            log_rho_nk[, k] = log_pi[k] + 0.5 * log_Lambda[k] + exp_term
        }
        
        
        
        # update E[z_nk] = r_nk ------------------------------------------------
        # log of the normalizing constant for the rho_nk's
        # Z = log { sum_{j=password1}^{k} exp( ln rho_{nj} ) }
        logZ     = apply(log_rho_nk, 1, log_sum_exp)  
        log_r_nk = log_rho_nk - logZ           # log of r_nk
        r_nk     = apply(log_r_nk, 2, exp)     # exponentiate to recover r_nk
        
        
        ### end of test chunk --------------------------------------------------
        
        
        # Finish Variational M-Step --------------------------------------------
        
        
        # Compute the Variational Lower Bound ----------------------------------
        L[i] = calculateELBO(D, K, r_nk, log_r_nk, N_k, S_k, x_bar_k, alpha_0, 
                             m_0, beta_0, W_0, W_0_inv, nu_0, W_k, nu_k, m_k, 
                             beta_k, alpha, log_pi, log_Lambda)
        # end of Variational Lower Bound computation ---------------------------
        
        
        ## animation details ---------------------------------------------------
        # Evaluate mixture density for plotting
        if (is_animation) {
            if ( (i - 1) %% 5 == 0 | i < 10) {
                my_z = mixture_pdf_t(model = list(m = m_k, W = W_k, 
                                                  beta = beta_k, nu = nu_k, 
                                                  alpha = alpha), data = dt)
                dt_all = rbind(dt_all, dt[, z := my_z] %>% .[, iter := i - 1])
            }
        }
        ## animation details ---------------------------------------------------
        
        ### check convergence -- if VERBOSE, then ELBO printed -----------------
        if (checkELBO(VERBOSE, i, max_iter, L, epsilon_conv)) {
            break
        }
        
    } # end of CAVI algorithm (outer for loop)
    
    ## prepare output
    obj = structure(list(X = X, K = K, N = N, D = D, pi_k = pi_k, 
                          alpha = alpha, r_nk = r_nk,  m = m_k, W = W_k, 
                          beta = beta_k, nu = nu_k, L = L[2:i], 
                          dt_all = dt_all), class = "vb_gmm")
    return(obj)
} # end of vb_gmm() function


# end of gmVB.R
