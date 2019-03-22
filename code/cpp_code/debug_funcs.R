
# debug_funcs.R


# checkCorrectness() : check equality (up to a certain tolerance) of 
#                      variational parameters calculated by R/C++ 
#                      implementation of e-step, m-step
# input : 
#        theta_R   : object w/ var. params calculated in R
#        theta_cpp : object w/ var. params calculated in C++
checkCorrectness = function(theta_R, theta_cpp) {
    
    # parmaeters whose types need to be changed so that they match:
    # Q_k, Q_k_inv, V_k, V_k_inv 
    
    # test explicit variational parameter equality
    if (!isTRUE(all.equal(theta_R$beta_k, theta_cpp$beta_k))) {
        print("beta_k's not equal")
    }
    if (!isTRUE(all.equal(theta_R$tau_k, theta_cpp$tau_k))) {
        print("tau_k's not equal")
    }
    if (!isTRUE(all.equal(theta_R$gamma_k, theta_cpp$gamma_k))) {
        print("tau_k's not equal")
    }
    
    # e-step updates
    if (!isTRUE(all.equal(theta_R$r_nk, theta_cpp$r_nk))) {
        print("r_nk's not equal")
    }
    if (!isTRUE(all.equal(theta_R$N_k, theta_cpp$N_k))) {
        print("N_k's not equal")
    }
    
    # m-step updates
    
    # bouchard parameters
    if (!isTRUE(all.equal(theta_R$alpha, theta_cpp$alpha))) {
        print("alpha's not equal")
    }
    if (!isTRUE(all.equal(theta_R$xi, theta_cpp$xi))) {
        print("xi's not equal")
    }
    if (!isTRUE(all.equal(theta_R$lambda, theta_cpp$lambda))) {
        print("lambda's not equal")
    }
    if (!isTRUE(all.equal(theta_R$phi, theta_cpp$phi))) {
        print("phi's not equal")
    }
    
    # q(gamma_k)
    if (!isTRUE(all.equal(theta_R$eta_k, theta_cpp$eta_k))) {
        print("eta_k's not equal")
    }
    if (!isTRUE(all.equal(theta_R$mu_k, theta_cpp$mu_k))) {
        print("mu_k's not equal")
    }
    
    # q(beta_k | tau_k)
    if (!isTRUE(all.equal(theta_R$zeta_k, theta_cpp$zeta_k))) {
        print("zeta_k's not equal")
    }
    if (!isTRUE(all.equal(theta_R$m_k, theta_cpp$m_k))) {
        print("m_k's not equal")
    }
    
    # q(tau_k)
    if (!isTRUE(all.equal(theta_R$a_k, theta_cpp$a_k))) {
        print("a_k's not equal")
    }
    if (!isTRUE(all.equal(theta_R$b_k, theta_cpp$b_k))) {
        print("b_k's not equal")
    }
    
    # check all precision matrices together: Q_k, Q_k_inv, V_k, V_k_inv
    for (k in 1:length(theta_R$Q_k)) {
        
        # ensure that the objects being compared are of same type
        # in R implementation, if D = 1, then each (1 x 1) matrix will be 
        # converted into a scalar, which cannot be compared to the (1 x 1)
        # matrix produced by C++ implementation 
        # (all.equal() requires type match)
        
        V_k_mat     = matrix(theta_R$V_k[,,k], nrow = nrow(theta_cpp$V_k[[k]]),
                                               ncol = ncol(theta_cpp$V_k[[k]]))
        V_k_inv_mat = matrix(theta_R$V_k[,,k], nrow = nrow(theta_cpp$V_k[[k]]),
                                               ncol = ncol(theta_cpp$V_k[[k]]))
        Q_k_mat     = matrix(theta_R$V_k[,,k], nrow = nrow(theta_cpp$Q_k[[k]]),
                                               ncol = ncol(theta_cpp$Q_k[[k]]))
        Q_k_inv_mat = matrix(theta_R$V_k[,,k], nrow = nrow(theta_cpp$Q_k[[k]]),
                                               ncol = ncol(theta_cpp$Q_k[[k]]))
        
        if (!isTRUE(all.equal(V_k_mat, theta_cpp$V_k[[k]]))) {
            print("V_k's not equal")
        }
        if (!isTRUE(all.equal(V_k_inv_mat, theta_cpp$V_k_inv[[k]]))) {
            print("V_k_inv's not equal")
        }
        if (!isTRUE(all.equal(Q_k_mat, theta_cpp$Q_k_mat[[k]]))) {
            print("Q_k's not equal")
        }
        if (!isTRUE(all.equal(Q_k_inv_mat, theta_cpp$Q_k_inv_mat[[k]]))) {
            print("Q_k_inv's not equal")
        }
    }
    
    
} # end of checkCorrectness() function


# end of debug_funcs.R
