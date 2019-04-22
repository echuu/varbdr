

setwd("C:/Users/chuu/varbdr/code")

source("globals.R")

setwd(HOME_DIR)
source(UNI_CDE) 
source(DENSITY)

# initialize parameters needed for varbdr_cpp()
intercept = TRUE
max_iter  = 1e4

N = 1e3  # number of observations
K = 2    # number of clusters
D = 2    # dimension of covaraites

synth_data_bimodal = r_binorm(N)                     # generate data
y      = synth_data_bimodal$y                        # vector of responses
X      = cbind(1, synth_data_bimodal$X)              # include intercept term
params = synth_data_bimodal$param_mat                # needed for plots

setwd("C:/Users/chuu/varbdr/code/cpp_code")
source("debug_funcs.R")                              # load debug functions
sourceCpp("getVarParams.cpp")                        # load C++ implementation
source(paste(COV_DEP, VARBDR, sep = '/'))            # load R implementation

# when initializating variational parameters, we randomly initialize the
# mean components for beta_k, gamma_k; we save these so that the R and C++ 
# implementations 'start' in the same place (for verification purposes)
meanParams = generateParams(y, X, N, D, K, intercept, max_iter)
(m_k  = meanParams$m_k)
(mu_k = meanParams$mu_k)

# run both cavi algorithms
theta_R  = varbdr(y = y, X = X, K, intercept = FALSE, m_k = m_k, mu_k = mu_k)    
theta_cpp = varbdr_cpp(y, X, N, D, K, intercept, max_iter)

# verify results of both cavi algs are same (ELBO will be slightly different)
checkCorrectness(theta_R, theta_cpp)

# overlay approximating densities


y_grid = seq(-5, 5, length.out = 1000)
overlayPlots = function(theta, K, den_label,
                        x = c(-0.6745, 0, 0.6745)) {
    
    n = length(theta)
    approx_d  = rep(list(py_bouch), n)
    p_list2 = xQuantileDensity(x, params, d_binorm,
                               approx_d, den_label, theta, K,
                               y_grid)
    
    p = multiplot(plotlist = as.list(sapply(p_list2, function(x) x[1])), 
                  cols = 3)
    return(p)
}


# intercept now only plays role in new input for plotting
# if intercept is fitted, first column of X should be vector of 1's 
# before passing into varbdr() and varbdr_cpp()
overlayPlots(theta = list(theta_cpp), K, den_label = c("N=1e4"))




# compare performance (takes a long time)
microbenchmark(varbdr_cpp(y, X, N, D, K, intercept, max_iter),
               varbdr(y = y, X = X, K, intercept, m_k = m_k, mu_k = mu_k),
               times = 20)




library("np")
X_np = data.frame(X[,2])
y_np = data.frame(y)
xy_df = data.frame(y = y_np, x = X_np)
names(xy_df) = c("y", "x")

fy.x = npcdens(y ~ x, xy_df)

x = c(-0.6745, 0, 0.6745)
uni_res = plotCD(theta_cpp, K, x, y_grid, d_binorm, k_den = fy.x)
uni_res$cd_plots










