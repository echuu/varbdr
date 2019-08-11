
# elbo.R --- implementation of the calculation of the variational lower bound



# computeELBO() : compute the variational lower bound using the current settings
#                 of the variational family, compare value of lower bound w/
#                 previous iteration's lower bound and determine convergence
#                 status
computeELBO = function(prior, theta) {
    
    
    
    
    
    # determine convergence status
    
    print(paste('iter:', 1, conv_status, sep = ' '))
    print(paste('iter:', 2, conv_status, sep = ' '))
    
    return(theta)
    
} # end computeELBO() function



# end of computeELBO.R 

