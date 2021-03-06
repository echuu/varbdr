

# test_modules.R

source("~/varbdr/code/globals.R")
setwd(HOME_DIR)
source(DP_BDR)    # so that we have access to data



setwd("~/varbdr/code/var_select")


setwd("~/varbdr/code/cov_indep")
source("misc.R")



# generate data
N = 500
D = 1
K = 2
synth_data_1d = r_dpmix2(N)
y = synth_data_1d$y
X = synth_data_1d$X
X = scale(X, scale = TRUE)

intercept = FALSE;     # dated setting? used to check that 1st column is 1 vec

# initialize parameters needed to pass into initPriors() function
# TODO: figure out the appropriate value to initialize pi_d with
# C&S implementation seems to normalize the probabilities over all features (?)
m_0 = 0; xi_0 = 0.3; pi_d = 0.7;                 # beta params             
a_0 = 1; b_0 = 1;                                # tau params
g_0 = 0; Sigma_0 = diag(rep(1, D));              # gamma params

VERBOSE = TRUE; tol = 1e-3;  max_iter = 1e4;     # algo/converge params


# test initPriors() function
source("initPriors.R")
prior = initPriors(y, X, K, 
                   m_0, xi_0, pi_d,              # beta components
                   a_0, b_0,                     # tau components
                   g_0, Sigma_0,                 # gamma components 
                   max_iter, tol, VERBOSE)

# test initVarParams() function
source("initVarParams.R")
theta = initVarParams(y, X, N, D, K, intercept = FALSE, max_iter, 
                      m_d = NULL, mu_k = NULL)

#### test updateFunctions.R
source("updateFunctions.R")
source("eStep_vs.R")
source("mStep_vs.R")



out_tau = precisionUpdate(prior, theta) # runs with no complaints

out_gamma = weightUpdate(prior, out_tau)

out_ss = spikeSlabUpdate(prior, out_gamma) # runs with no complaints




source("updateFunctions.R")
theta_eStep = eStep_vs(prior, theta)

out_tau = precisionUpdate(prior, theta_eStep) # runs with no complaints

out_gamma = weightUpdate(prior, out_tau)

out_ss = spikeSlabUpdate(prior, out_gamma) # runs with no complaints



theta_mStep = mStep_vs(prior, theta_eStep)


# will return NaN for now because r_nk = 0; need to implement e-step before this
# (hopefully) returns real values
source("computeELBO.R")
elbo_0 = elbo_vs(prior, theta_mStep)

elbo_0
                                        

