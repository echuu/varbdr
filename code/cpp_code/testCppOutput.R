

# source("C:/Users/chuu/varbdr/code/globals.R")
source("~/varbdr/code/globals.R")
setwd(HOME_DIR)
source(DP_BDR) 
source(DENSITY)

# setwd("C:/Users/chuu/varbdr/code/cpp_code")
setwd("~/varbdr/code/cpp_code")

# source("C:/Users/chuu/varbdr/code/cov_dep/misc.R")
source("~/varbdr/code/cov_dep/misc.R")
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
theta = initVarParams(y, X, N, D, K, intercept, max_iter)

theta = mStep(theta, prior)

# compare corresponding componentes in theta_cpp and theta

sourceCpp("getVarParams.cpp")
theta_cpp = testConstructor(y, X, N, D, K, intercept, max_iter)
str(theta_cpp)

theta_cpp$alpha
theta_cpp$xi
theta_cpp$phi
theta_cpp$lambda

theta$alpha
theta$xi
theta$phi
theta$lambda

theta_cpp$Q_k
theta$Q_k

all.equal(theta_cpp$alpha, theta$alpha)
all.equal(theta$xi, theta_cpp$xi)
all.equal(theta$phi, theta_cpp$phi)
all.equal(theta$lambda, theta_cpp$lambda)

# 2nd half of m-step testing
checkEqual(theta_cpp$Q_k_inv, theta$Q_k_inv)
all.equal(theta_cpp$V_k_inv, theta$V_k_inv)

all.equal(theta_cpp$m_k, theta$m_k)
all.equal(theta_cpp$mu_k, theta$mu_k)


# TODO: test these after incorporating the rest of the
#       m-step code

all.equal(theta_cpp$a_k, theta$a_k)
all.equal(theta_cpp$b_k, theta$b_k)


# one iteration of m-step performance
microbenchmark(testConstructor(y, X, N, D, K, intercept, max_iter), 
               mStep(theta, prior))








