

source("C:/Users/chuu/varbdr/code/globals.R")
# source("~/varbdr/code/globals.R")
setwd(HOME_DIR)
source(DP_BDR) 
source(DENSITY)

setwd("C:/Users/chuu/varbdr/code/cpp_code")
# setwd("~/varbdr/code/cpp_code")

source("C:/Users/chuu/varbdr/code/cov_dep/misc.R")
# source("~/varbdr/code/cov_dep/misc.R")
source("debug_funcs.R")


N = 100
K = 3
D = 1

synth_data_1d = r_dpmix2(N)

y      = synth_data_1d$y
X      = synth_data_1d$X
intercept = FALSE
max_iter = 5000

sourceCpp("getVarParams.cpp")
theta_cpp = testConstructor(y, X, N, D, K, intercept, max_iter)
str(theta_cpp)

mu_k = theta_cpp$mu_k
m_k  = theta_cpp$m_k

# model parameters to pass into R implementation of varbdr() -------------------
# note: these are defined internally in C++ implementation

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

# use c++ initializations for these variational parameters


prior = initPriors(y, X, K, m_0, Lambda_0, a_0, b_0, g_0, Sigma_0,
                   max_iter, tol, VERBOSE)

# set back to initVarParams_0() later
theta = initVarParams_0(y, X, N, D, K, intercept, max_iter, m_k, mu_k)

# check equality after initialization
checkCorrectness(theta, theta_cpp)

theta = eStep(theta, prior)
theta = mStep(theta, prior)
str(theta)
checkCorrectness(theta, theta_cpp)


# compare corresponding componentes in theta_cpp and theta

sourceCpp("getVarParams.cpp")
theta_cpp = testConstructor(y, X, N, D, K, intercept, max_iter)
str(theta_cpp)



# one complete iteration of cavi: eStep(), mStep()

source(paste(COV_DEP,  E_STEP,          sep = '/'))
cavi_R = function() {
    prior = initPriors(y, X, K, m_0, Lambda_0, a_0, b_0, g_0, Sigma_0,
                       max_iter, tol, VERBOSE)
    theta = initVarParams_0(y, X, N, D, K, intercept, max_iter, m_k, mu_k)
    for (i in 1:20) {
        theta = eStep(theta, prior)
        theta = mStep(theta, prior)
        theta = elbo(theta, prior)
    }
}

sourceCpp("getVarParams.cpp")
cavi_cpp = function() {
    theta_cpp = testConstructor(y, X, N, D, K, intercept, max_iter)
}


### test elbo
source(paste(COV_DEP,  ELBO,            sep = '/'))
theta = initVarParams_0(y, X, N, D, K, intercept, max_iter, m_k, mu_k)
theta = eStep(theta, prior)
theta = mStep(theta, prior)
theta = elbo(theta, prior)

theta$e_ln_p_y           # --- match
theta$e_ln_p_z           # --- match
theta$e_ln_p_gamma       # --- match
theta$e_ln_p_beta_tau    # --- match
theta$e_ln_q_z           # --- match
theta$e_ln_q_beta_tau    # --- match
theta$e_ln_q_gamma       # --- match


sourceCpp("getVarParams.cpp")
theta_cpp = testConstructor(y, X, N, D, K, intercept, max_iter)





