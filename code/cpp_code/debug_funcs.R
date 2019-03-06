
# debug_funcs.R


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
