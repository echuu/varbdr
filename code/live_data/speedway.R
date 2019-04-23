
setwd("C:/Users/chuu/varbdr/code")

source("globals.R")

setwd(HOME_DIR)
source(UNI_CDE) 
source(DENSITY)
source(paste(COV_DEP, VARBDR, sep = '/'))            # load R implementation


setwd("C:/Users/chuu/varbdr/code/cpp_code")
source("debug_funcs.R")                              # load debug functions
sourceCpp("getVarParams.cpp")                        # load C++ implementation
source(paste(COV_DEP, VARBDR, sep = '/'))            # load R implementation


library(hdrcde)
library(np)

# 
# lane2 : (1318 x 2) lane 2 of the 4-lane California freeway I-880
# flow  : traffic flow in vehicles per lane per hour
# speed : speed in miles per hour

# plot of traffic speed vs. traffic flow
plot(lane2) 
lane2_sub = lane2 %>% filter(flow >= 1000, flow <= 1620)
plot(lane2_sub, main = "1000 < flow < 1620")
abline(v = 1400, col = 'red')

# format data
X = cbind(as.matrix(lane2_sub$flow))
X = X / 100
X = cbind(1, X)      # fit intercept
y = lane2_sub$speed
N = nrow(X)
D = ncol(X)
K = 2                # bimodal by observation
intercept = TRUE
max_iter  = 1e5

x_in = c(1400 / 100)

theta_sf = varbdr_cpp(y, X, N, D, K, intercept, max_iter) # converge in 4809

# meanParams = generateParams(y, X, N, D, K, intercept = FALSE, max_iter)
# (m_k  = meanParams$m_k)
# (mu_k = meanParams$mu_k)
# theta_R  = varbdr(y = y, X = X, K, intercept = FALSE, m_k = m_k, mu_k = mu_k)    
# N_evals    = length(y_grid)
# data_ygrid = list(y = y_grid, x = c(1, x_in / 100))

# obtain kernel-based estimate
X_np = data.frame(X[,2])
y_np = data.frame(y)
xy_df = data.frame(y = y_np, x = X_np)
names(xy_df) = c("y", "x")

fy.x = npcdens(y ~ x, xy_df)

# plot entire density
y_grid = seq(min(y), max(y) + 20, length.out = 2000)
sf_res = plotCD(theta_sf, K, x_in, y_grid, true_den = NULL, k_den = fy.x)
sf_res$cd_plots[[1]] + ggtitle("vb (red), kernel (blue)")

# examine the location of where the 2nd mode should be
y_grid = seq(20, 50, length.out = 2000)
sf_res = plotCD(theta_sf, K, x_in/100, y_grid, true_den = NULL, k_den = fy.x)
sf_res$cd_plots[[1]] + ggtitle("vb (red), kernel (blue)")


# ------------------------------------------------------------------------------

#### ISLR data

library(ISLR)

table(Wage$year) # 2003-2009
table(Wage$age)


plot(Wage$age,Wage$logwage)

plot(density(Wage$logwage[Wage$age == 34]))
plot(density(Wage$logwage[Wage$age == 42]))
plot(density(Wage$logwage[Wage$age == 51]))

h <- hist(Wage$logwage[Wage$age == 34], breaks = 9, plot=FALSE)
h$counts=h$counts/sum(h$counts)
plot(h)

par(mfrow = c(1,3))
hist(Wage$logwage[Wage$age == 34], main = "age = 34")
hist(Wage$logwage[Wage$age == 42], main = "age = 42")
hist(Wage$logwage[Wage$age == 51], main = "age = 51")
par(mfrow = c(1, 1))

x_in = c(34, 42, 51) / 10

# X = as.matrix(as.double(Wage$age[(Wage$age >= 32) & (Wage$age <= 58)])) / 10
# y = Wage$logwage[(Wage$age >= 32) & (Wage$age <= 58)]

X = as.matrix(as.double(Wage$age)) / 10
y = Wage$logwage


N = nrow(X)
D = ncol(X)
K = 2
intercept = FALSE
max_iter  = 7e3

y_grid = seq(min(y) - 2.5, max(y) + 3, length.out = 2000)

theta_wage = varbdr_cpp(y, X, N, D, K, intercept, max_iter) # 1160

# compute density estimate
X_np = data.frame(X[,1])
y_np = data.frame(y)
xy_df = data.frame(y = y_np, x = X_np)
names(xy_df) = c("y", "x")

fyx_wage = npcdens(y ~ x, xy_df)

wage_cd = plotCD(theta_wage, K, x_in, y_grid, true_den = NULL, k_den = NULL)

multiplot(plotlist = wage_cd$cd_plots, cols = 3)




















###### debugging

m_0 = c(colMeans(X))                         # normal params              
Lambda_0 = diag(rep(1, ncol(X))) 
a_0 = 1
b_0 = 1                             # gamma params
g_0 = 0 
Sigma_0 = diag(rep(1, ncol(X)))


###### debugging

source(paste(COV_DEP, VARBDR, sep = '/'))            # load R implementation

prior = initPriors(y, X, K, m_0, Lambda_0, a_0, b_0, g_0, Sigma_0,
                   max_iter, tol, VERBOSE)

theta = initVarParams(y, X, N, D, K, intercept, max_iter, m_k, mu_k)

theta = eStep(theta, prior)
theta = mStep(theta, prior)
theta = elbo(theta, prior)

# check xi (previously had too large values), phi (overflow in exponential term)
head(theta$xi, 50)
head(theta$phi, 50)

# check the expectations
theta$e_ln_p_y
theta$e_ln_p_z
theta$e_ln_p_gamma
theta$e_ln_p_beta_tau
theta$e_ln_q_z
theta$e_ln_q_beta_tau
theta$e_ln_q_gamma










