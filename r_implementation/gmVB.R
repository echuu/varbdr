library(matrixcalc)
library(ggplot2)
library(data.table)
library(purrr)        # overlapping transpose function from data.table
library(mvtnorm)
library(Matrix)


setwd("~/varbdr/r_implementation")

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
    # W_0     : prior scale matrix for lambda
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
    # W       : updated scale matrix for lambda
    # nu      : updated degrees of freedom for lambda
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
    r_nk      = log_r_nk = log_rho_nk = matrix(0, nrow = N, ncol = K)
    x_bar_k   = matrix(0, nrow = D, ncol = K)         # Bishop 10.52
    S_k       = W_k = array(0, c(D, D, K ))           # Bishop 10.53
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
        
        ##-------------------------------
        # Variational E-Step
        ##-------------------------------
        
        # update responsibilities: involves E[log(det(Lambda_k))], 
        # E[log(det(pi_k))], E[(x_n - mu_k)' W_k (x_n - mu_k)]
        for (k in 1:K) {
            # vectorize the operation: (x_n - mu_k), n = 1,...,N
            # (N X D) matrix differences stored row-wise
            diff = sweep(X, MARGIN = 2, STATS = m_k[, k], FUN = "-")
            exp_term = - 0.5 * D / beta_k[k] - 0.5 * nu_k[k] * 
                diag(diff %*% W_k[,,k] %*% t(diff))
            
            # log of 10.67 (done for k = 1,...,K simultaneously)
            log_rho_nk[, k] = log_pi[k] + 0.5 * log_Lambda[k] + exp_term
        }
        
        # Responsibilities using the logSumExp for numerical stability
        Z        = apply(log_rho_nk, 1, log_sum_exp)
        log_r_nk = log_rho_nk - Z              # log of 10.49
        r_nk     = apply(log_r_nk, 2, exp)     # 10.49
        
        ##-------------------------------
        # Variational M-Step
        ##-------------------------------
        N_k = colSums(r_nk) + 1e-10  # K x 1 : N_1, N_2, ... , N_K ----  10.51
        for (k in 1:K) {
            x_bar_k[, k] = (r_nk[ ,k] %*% X) / N_k[k]   # 10.52
            x_cen        = sweep(X,MARGIN = 2,STATS = x_bar_k[, k],FUN = "-")
            S_k[, , k]   = t(x_cen) %*% (x_cen * r_nk[, k]) / N_k[k]  # 10.53
        }
        
        ## Update Dirichlet parameter
        alpha = alpha_0 + N_k  # 10.58
        ## Compute expected value of mixing proportions: E(pi_k)
        pi_k = (alpha_0 + N_k) / (K * alpha_0 + N)
        
        ## Update parameters for Gaussian-Wishart distribution
        beta_k = beta_0 + N_k      # 10.60
        nu_k   = nu_0 + N_k + 1    # 10.63
        for (k in 1:K) {
            # 10.61
            m_k[, k] = (1 / beta_k[k]) * (beta_0 * m_0 + N_k[k] * x_bar_k[, k])  
            # 10.62
            W_k[, , k] = W_0_inv + N_k[k] * S_k[,,k] + 
                ((beta_0 * N_k[k]) / (beta_0 + N_k[k])) * 
                tcrossprod((x_bar_k[, k] - m_0))    
            W_k[, , k] = solve(W_k[, , k])
        }
        ## Update expectations over \pi and \Lambda
        # 10.66
        log_pi = digamma(alpha) - digamma(sum(alpha))                      
        for (k in 1:K) { # 10.65                                              
            log_Lambda[k] = sum(digamma((nu_k[k] + 1 - 1:D) / 2)) + 
                D * log(2) + log(det(W_k[,,k])) 
        }
        
        ##-------------------------------
        # Variational lower bound
        ##-------------------------------
        lb_px = lb_pml = lb_pml2 = lb_qml = 0
        for (k in 1:K) {
            
            # 10.71
            lb_px = lb_px + N_k[k] * 
                (log_Lambda[k] - D / beta_k[k] - nu_k[k] * 
                     matrix.trace(S_k[,,k] %*% W_k[,,k]) - 
                     nu_k[k] * t(x_bar_k[, k] -  m_k[,k]) %*% W_k[,,k] %*% 
                     (x_bar_k[, k] - m_k[, k]) - D * log(2 * pi) ) 
            
            # 10.74
            lb_pml = lb_pml + D * log(beta_0  /(2 * pi)) + log_Lambda[k] - 
                (D * beta_0) / beta_k[k] - beta_0 * nu_k[k] * 
                t(m_k[,k] - m_0) %*% W_k[,,k] %*% (m_k[,k] - m_0)    
            
            # 10.74
            lb_pml2 = lb_pml2 + nu_k[k] * matrix.trace(W_0_inv %*% W_k[,,k]) 
            
            # 10.77
            lb_qml = lb_qml + 0.5 * log_Lambda[k] + 
                0.5 * D * log(beta_k[k] / (2 * pi)) - 
                0.5 * D - logB(W = W_k[,,k], nu = nu_k[k]) - 
                0.5 * (nu_k[k] - D - 1) * log_Lambda[k] + 0.5 * nu_k[k] * D
        }
        
        
        # ELBO computation --- move this to another file -----------------------
        lb_px  = 0.5 * lb_px                                            # 10.71
        lb_pml = 0.5 * lb_pml + K * logB(W = W_0, nu = nu_0) + 0.5 * 
            (nu_0 - D - 1) * sum(log_Lambda) - 0.5 * lb_pml2            # 10.74
        lb_pz  = sum(r_nk %*% log_pi)                                   # 10.72
        lb_qz  = sum(r_nk * log_r_nk)                                   # 10.75
        lb_pp  = sum((alpha_0 - 1) * log_pi) + lgamma(sum(K * alpha_0)) -
            K * sum(lgamma(alpha_0))                                    # 10.73
        lb_qp  = sum((alpha - 1) * log_pi) + lgamma(sum(alpha)) - 
            sum(lgamma(alpha))                                          # 10.76
        
        # Sum all parts to compute lower bound
        L[i] = lb_px + lb_pz + lb_pp + lb_pml - lb_qz - lb_qp - lb_qml
        
        # end of ELBO computation ----------------------------------------------
        
        
        
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
        
        ### check convergence -- if VERBOSE, then ELBO printed
        if (checkELBO(VERBOSE, i, max_iter, L, epsilon_conv)) {
            break
        }
        
    } # end of CAVI algorithm
    
    ## prepare output
    obj = structure(list(X = X, K = K, N = N, D = D, pi_k = pi_k, 
                          alpha = alpha, r_nk = r_nk,  m = m_k, W = W_k, 
                          beta = beta_k, nu = nu_k, L = L[2:i], 
                          dt_all = dt_all), class = "vb_gmm")
    return(obj)
} # end of vb_gmm() function


# end of gmVB.R
