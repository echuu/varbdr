library(matrixcalc)
library(ggplot2)
library(data.table)
library(purrr)        # overlapping transpose function from data.table
library(mvtnorm)
library(Matrix)


# setwd("~/varbdr/r_implementation")
setwd("C:/Users/chuu/varbdr/cov_independent_code/r_implementation")

source("densityCalculations.R")
source("vb_calcs.R")
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
    
    # Iterate to find optimal parameters
    for (i in 2:max_iter) {
        
        ## Variational E-Step
        theta = eStep(N, D, K, X, theta)
        
        ## Variational M-Step
        theta = mStep(N, K, D, X, theta, alpha_0, beta_0, nu_0, W_0_inv, m0)
        
        
        # Compute the Variational Lower Bound ----------------------------------
        L[i] = elbo(D, K, theta, alpha_0, m_0, beta_0, W_0, W_0_inv, nu_0)
        # end of Variational Lower Bound computation ---------------------------
        
        
        ## animation details ---------------------------------------------------
        # Evaluate mixture density for plotting
        if (is_animation) {
            if ( (i - 1) %% 5 == 0 | i < 10) {
                my_z = mixture_pdf_t(model = list(m = theta$m_k, W = theta$W_k, 
                                                  beta = theta$beta_k, 
                                                  nu = theta$nu_k, 
                                                  alpha = theta$alpha),
                                     data = dt)
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
    obj = structure(list(X = X, K = K, N = N, D = D, pi_k = theta$pi_k, 
                         alpha = theta$alpha, r_nk = theta$r_nk,  m = theta$m_k, 
                         W = theta$W_k, beta = theta$beta_k, nu = theta$nu_k, 
                         L = L[2:i], dt_all = dt_all), class = "vb_gmm")
    return(obj)
} # end of vb_gmm() function


# end of gmVB.R
