

setwd("C:/Users/chuu/varbdr/code")

source("globals.R")

setwd(HOME_DIR)
source(UNI_CDE) 
source(DENSITY)


N = 500
K = 3
D = 1

synth_data_bimodal = r_binorm(N)
y      = synth_data_bimodal$y
X      = synth_data_bimodal$X
params = synth_data_bimodal$param_mat
y_grid = seq(-5, 5, length.out = 1000)

x         = c(-0.6745, 0, 0.6745)
f_yx      = matrix(0, nrow = length(y_grid), ncol = length(x))

for(i in 1:length(x)) {
    f_yx[,i] = d_binorm(y_grid = y_grid, x = x[i])
}

df_y      = data.frame(x = y_grid, f_yx)
df_y_long = melt(df_y, id.vars = "x")
levels(df_y_long$variable) = c("x = -0.6745", "x = 0", "x = 0.6745")
ggplot(df_y_long, aes(x, y = value)) + geom_point(size = 0.9) +
    scale_x_continuous(breaks = seq(min(y_grid), max(y_grid), by = 1)) + 
    facet_wrap(~variable)

X = cbind(1,X)
D = 2
# run cavi ---------------------------------------------------------------------

intercept = FALSE
max_iter = 1000


setwd("C:/Users/chuu/varbdr/code/cpp_code")
source("debug_funcs.R")
sourceCpp("getVarParams.cpp")
meanParams = generateParams(y, X, N, D, K, intercept, max_iter)
(m_k  = meanParams$m_k)
(mu_k = meanParams$mu_k)

source(paste(COV_DEP,  INIT_PRIORS,     sep = '/'))
source(paste(COV_DEP,  INIT_VAR_PARAMS, sep = '/'))
source(paste(COV_DEP,  E_STEP,          sep = '/'))
source(paste(COV_DEP,  M_STEP,          sep = '/'))
source(paste(COV_DEP,  ELBO,            sep = '/'))
source(paste(COV_DEP,  MISC_FUNCS,      sep = '/'))
source(paste(HOME_DIR, DENSITY,         sep = '/'))

m_0 = c(colMeans(X))                         # normal params              
Lambda_0 = diag(rep(1, ncol(X))) 
a_0 = 1
b_0 = 1                             # gamma params
g_0 = 0 
Sigma_0 = diag(rep(1, ncol(X)))
tol = 1e-3
VERBOSE = FALSE

prior = initPriors(y, X, K, m_0, Lambda_0, a_0, b_0, g_0, Sigma_0,
                   max_iter, tol, VERBOSE)
theta = initVarParams_0(y, X, N, D, K, intercept, max_iter, 
                        m_k = m_k, mu_k = mu_k)
checkCorrectness(theta, meanParams)                      # correct up to here

theta = eStep(theta, prior)                              # correct up to here
theta = mStep(theta, prior)                              # correct up to here
theta = elbo(theta, prior)   # checking

sourceCpp("getVarParams.cpp")
theta1_cpp = varbdr_cpp(y, X, N, D, K, intercept, max_iter) # 1 iter estep
checkCorrectness(theta, theta1_cpp)


# ------------------------------------------------------------------------------



sourceCpp("getVarParams.cpp")
source("debug_funcs.R")
checkCorrectness(theta, meanParams)

theta1_cpp = varbdr_cpp(y, X_cpp, N, D, K, intercept, max_iter)
checkCorrectness(theta, theta_cpp)




# ------------------------------------------------------------------------------



# run R implementation of CAVI
source(paste(COV_DEP, VARBDR, sep = '/'))            # load cov-dep varbdr.R
theta1_R = varbdr(y = y, X = X, K, intercept = FALSE, m_k = m_k, mu_k = mu_k)    

# run C++ implementation of CAVI
sourceCpp("getVarParams.cpp")
theta1_cpp = varbdr_cpp(y, X, N, D, K, intercept, max_iter)


source("debug_funcs.R")
checkCorrectness(theta1_R, theta1_cpp)


# overlay approximation and true density
overlayPlots = function(theta, K, den_label,
                        x = c(-0.6745, 0, 0.6745)) {
    
    # source(DENSITY)
    
    n = length(theta)
    approx_d  = rep(list(py_bouch), n)
    # approx_d = list(py_bouch, py_bouch, py_bouch, py_bouch)
    # den_label = c("cov-dep-1", "cov-dep-2", "cov-dep-3", "cov-dep-4")
    p_list2 = xQuantileDensity(x, params, d_binorm,
                               approx_d, den_label, theta, K,
                               y_grid)
    
    p = multiplot(plotlist = as.list(sapply(p_list2, function(x) x[1])), 
                  cols = 3)
    return(p)
}

overlayPlots(theta = list(theta1_R), 
             K = 2, den_label = c("N=500"))


