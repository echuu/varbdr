
source("C:/Users/chuu/varbdr/code/globals.R")

setwd(HOME_DIR)
source(DP_BDR) 
source(DENSITY)
source(paste(COV_DEP, VARBDR, sep = '/')) # load cov-dep varbdr.R

setwd("C:/Users/chuu/varbdr/code/cpp_code")
source("debug_funcs.R")

N = 500
K = 2
synth_data_1d = r_dpmix2(N)
y = synth_data_1d$y
X = synth_data_1d$X
y_grid = seq(-3, 1.5, length.out = 1000)

sourceCpp("getVarParams.cpp")
meanParams = generateParams(y, X, N, D, K, intercept, max_iter)
m_k  = meanParams$m_k
mu_k = meanParams$mu_k 


# run both cavi algorithms
theta0_R = varbdr(y, X, K, m_k = m_k, mu_k = mu_k)   # store CAVI results

theta0_cpp = varbdr_cpp(y, X, N, D, K, intercept, max_iter)

# check equality after convergence
checkCorrectness(theta0_R, theta0_cpp)

# compare speed
microbenchmark(varbdr(y, X, K, m_k = m_k, mu_k = mu_k),
               varbdr_cpp(y, X, N, D, K, intercept, max_iter),
               times = 25)




## check plots -- need to re-edit the functions to take C++ return

theta = list(theta0_R)
overlayPlots = function(theta, K, den_label,
                        x = c(0.15, 0.25, 0.49, 0.75, 0.88, 0.95)) {
    
    # source(DENSITY)
    
    n = length(theta)
    approx_d  = rep(list(py_bouch), n)
    p_list2 = xQuantileDensity(x, params, d_dpmix2,
                               approx_d, den_label, theta, K,
                               y_grid)
    
    p = multiplot(plotlist = as.list(sapply(p_list2, function(x) x[1])), 
                  cols = 3)
    return(p)
}


overlayPlots(theta, K = 2, 
             den_label = c("N=500"))

