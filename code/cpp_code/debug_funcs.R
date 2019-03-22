
# debug_funcs.R


# checkCorrectness() : check equality (up to a certain tolerance) of 
#                      variational parameters calculated by R/C++ 
#                      implementation of e-step, m-step
# input : 
#        theta_R   : object w/ var. params calculated in R
#        theta_cpp : object w/ var. params calculated in C++
checkCorrectness = function(theta_R, theta_cpp) {
    
    
    
    
}



checkEqual = function(obj1, obj2) {
    VALID = TRUE
    
    for (i in 1:length(obj1)) {
        if (!all.equal(obj1[[i]], obj2[[i]])) {
            print(paste("object ", i, " does not match", sep = ''))
            VALID = FALSE
        }
    
        if (is.null(obj1[[1]]) || is.null(obj2[[2]])) {
            print(paste("object ", i, " is NULL", sep = ''))
            VALID = FALSE
        }
    }
    
    if (VALID) {
        print("Both objects equal.")
    }
    
} # end checkEqual() function





# end of debug_funcs.R
