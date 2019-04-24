
# debug_funcs.R


# checkCorrectness() : check equality (up to a certain tolerance) of 
#                      variational parameters calculated by R/C++ 
#                      implementation of e-step, m-step
# input : 
#        theta_R   : object w/ var. params calculated in R
#        theta_cpp : object w/ var. params calculated in C++
checkCorrectness = function(theta_R, theta_cpp) {
    
    
    ALL_CORRECT = TRUE
    
    # parmaeters whose types need to be changed so that they match:
    # Q_k, Q_k_inv, V_k, V_k_inv 
    
    # test explicit variational parameter equality
    if (!isTRUE(all.equal(theta_R$beta_k, theta_cpp$beta_k))) {
        ALL_CORRECT = FALSE
        print("beta_k's not equal")
    }
    if (!isTRUE(all.equal(theta_R$tau_k, theta_cpp$tau_k))) {
        ALL_CORRECT = FALSE
        print("tau_k's not equal")
    }
    if (!isTRUE(all.equal(theta_R$gamma_k, theta_cpp$gamma_k))) {
        ALL_CORRECT = FALSE
        print("tau_k's not equal")
    }
    
    # e-step updates
    if (!isTRUE(all.equal(theta_R$r_nk, theta_cpp$r_nk))) {
        ALL_CORRECT = FALSE
        print("r_nk's not equal")
    }
    # if (!isTRUE(all.equal(theta_R$N_k, theta_cpp$N_k))) {
    #     print("N_k's not equal")
    # }
    
    # m-step updates
    
    # bouchard parameters
    if (!isTRUE(all.equal(theta_R$alpha, theta_cpp$alpha))) {
        ALL_CORRECT = FALSE
        print("alpha's not equal")
    }
    if (!isTRUE(all.equal(theta_R$xi, theta_cpp$xi))) {
        ALL_CORRECT = FALSE
        print("xi's not equal")
    }
    if (!isTRUE(all.equal(theta_R$lambda, theta_cpp$lambda))) {
        ALL_CORRECT = FALSE
        print("lambda's not equal")
    }
    if (!isTRUE(all.equal(theta_R$phi, theta_cpp$phi))) {
        ALL_CORRECT = FALSE
        print("phi's not equal")
    }
    
    # q(gamma_k)
    if (!isTRUE(all.equal(theta_R$eta_k, theta_cpp$eta_k))) {
        ALL_CORRECT = FALSE
        print("eta_k's not equal")
    }
    if (!isTRUE(all.equal(theta_R$mu_k, theta_cpp$mu_k))) {
        ALL_CORRECT = FALSE
        print("mu_k's not equal")
    }
    
    # q(beta_k | tau_k)
    if (!isTRUE(all.equal(theta_R$zeta_k, theta_cpp$zeta_k))) {
        ALL_CORRECT = FALSE
        print("zeta_k's not equal")
    }
    if (!isTRUE(all.equal(theta_R$m_k, theta_cpp$m_k))) {
        ALL_CORRECT = FALSE
        print("m_k's not equal")
    }
    
    # q(tau_k)
    if (!isTRUE(all.equal(theta_R$a_k, theta_cpp$a_k))) {
        ALL_CORRECT = FALSE
        print("a_k's not equal")
    }
    if (!isTRUE(all.equal(theta_R$b_k, theta_cpp$b_k))) {
        ALL_CORRECT = FALSE
        print("b_k's not equal")
    }
    
    # check all precision matrices together: Q_k, Q_k_inv, V_k, V_k_inv
    for (k in 1:dim(theta_R$Q_k)[1]) {
        
        # ensure that the objects being compared are of same type
        # in R implementation, if D = 1, then each (1 x 1) matrix will be 
        # converted into a scalar, which cannot be compared to the (1 x 1)
        # matrix produced by C++ implementation 
        # (all.equal() requires type match)
        
        Vk     = matrix(theta_R$V_k[,,k], nrow = nrow(theta_cpp$V_k[[k]]),
                                          ncol = ncol(theta_cpp$V_k[[k]]))
        Vk_inv = matrix(theta_R$V_k_inv[,,k], nrow = nrow(theta_cpp$V_k[[k]]),
                                              ncol = ncol(theta_cpp$V_k[[k]]))
        Qk     = matrix(theta_R$Q_k[,,k], nrow = nrow(theta_cpp$Q_k[[k]]),
                                          ncol = ncol(theta_cpp$Q_k[[k]]))
        Qk_inv = matrix(theta_R$Q_k_inv[,,k], nrow = nrow(theta_cpp$Q_k[[k]]),
                                              ncol = ncol(theta_cpp$Q_k[[k]]))
        
        if (!isTRUE(all.equal(Vk, theta_cpp$V_k[[k]]))) {
            ALL_CORRECT = FALSE
            print("V_k's not equal")
        }
        if (!isTRUE(all.equal(Vk_inv, theta_cpp$V_k_inv[[k]]))) {
            ALL_CORRECT = FALSE
            print("V_k_inv's not equal")
        }
        if (!isTRUE(all.equal(Qk, theta_cpp$Q_k[[k]]))) {
            ALL_CORRECT = FALSE
            print("Q_k's not equal")
        }
        if (!isTRUE(all.equal(Qk_inv, theta_cpp$Q_k_inv[[k]]))) {
            ALL_CORRECT = FALSE
            print("Q_k_inv's not equal")
        }
    }
    
    if (ALL_CORRECT) {
        print("All tests passed! You a thug.")
    }
    
    
} # end of checkCorrectness() function


# end of debug_funcs.R
