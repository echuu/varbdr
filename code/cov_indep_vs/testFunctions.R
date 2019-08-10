

source("~/varbdr/code/globals.R")
setwd(HOME_DIR)
source(DP_BDR)    # so that we have access to data

setwd("~/varbdr/code/cov_indep_vs")
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

alpha_0 = rep(1 / K, K);                         # pi_k params
m_0 = 0; xi_0 = 0.3;                             # beta params             
pi_d = 0.3;                                      # prior prob of inclusion
a_0 = 1; b_0 = 1;                                # tau params
VERBOSE = TRUE; tol = 1e-3;  max_iter = 1e4;     # algo/converge params


# test initPriors() function
source("initPr.R")
prior = initPriors(y, X, K, 
                   alpha_0,                    # pi   : concentration param
                   m_0, xi_0,                  # beta : mean, scale for slab
                   pi_d,                       # beta : prior pip
                   a_0, b_0,                   # tau  : shape, rate
                   tol, VERBOSE)

# test initVarParams() function
source("initVP.R")
theta = initVarParams_indep(y, X, N, D, K)

source('updateFunctions.R')

# (0) test rnkUpdate() function
theta = rnkUpdate(prior, theta)

# (1) test ssUpdate() function
theta = ssUpdate(prior, theta)


# (2) test precUpdate() function
theta = precUpdate(prior, theta)

# (3) test wtUpdate() function
theta = wtUpdate(prior, theta)








