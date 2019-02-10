

source("globals.R")

setwd(HOME_DIR)
source(DP_BDR) 
source(DENSITY)


N = 500
K = 4
synth_data_1d = r_dpmix2(N)

y      = synth_data_1d$y
X      = synth_data_1d$X
params = synth_data_1d$param_mat

# 6th sublot in Figure 3
df_yx = data.frame(x = X[,2], y = y)
ggplot(df_yx, aes(x = x, y = y)) + geom_point() + stat_smooth() + theme_bw()


# 1st subplot
y_grid = seq(-0.5, 1.5, length.out = 500)
x = 0.49
f_yx = d_dpmix2(y_grid, c(1, x), c(x, x^4))
df_y = data.frame(x = y_grid, y = f_yx)
ggplot(df_y, aes(x, y)) + geom_point(size = 0.9) 

# run cavi ---------------------------------------------------------------------

# (1.1) covariate-INDEPENENT vb
source(paste(COV_INDEP, VARBDR, sep = '/'))      # load cov-indep varbdr.R
theta1_1 = varbdr(y = y, X = X, K)               # store CAVI results

# (1.2) covariate-DEPENDENT vb
source(paste(COV_DEP, VARBDR, sep = '/'))        # load cov-dep varbdr.R
theta1_2 = varbdr(y = y, X = X, K)               # store CAVI results


# generate plots ---------------------------------------------------------------

theta1      = list(theta1_1, theta1_2)      # list of var. params for each alg
DATA_GEN_ID = DP_MIX2

# generate plot overlays of the percentiles
source(DENSITY)
x = c(0.10, 0.25, 0.49, 0.75, 0.88)
approx_d  = list(py_0, py_bouch)            # list of approx density functions
den_label = c("cov-indep", "cov-dep (b)")   # labels for each approx density
p_list1 = xQuantileDensity(x, params, d_dpmix2,
                           approx_d, den_label, theta1, K,
                           DP_MIX2,
                           y_grid)

p_list1[[3]]


multiplot(plotlist = p_list1, cols = 2)


# first simutlation ------------------------------------------------------------


# fig 2 in DP : predictive density y_{n+1} at the 10, 25, 50, 75, 90 percentiles
#               of the empirical distribution of x_{i2}

# plot f(y|x), the true densities for varying values of x (0.1, 0.25, 0.5, 0.75)
#              to gain understanding of how the conditional density varies
#              as x goes from low to high

dp1_synth = r_dpmix1(N)
y0 = dp1_synth$y
X0 = dp1_synth$X
mu1 = dp1_synth$param_vec
y_grid = seq(-1.5, 1.5, length.out = 500)
x = 0.90

f_yx = d_dpmix1(y_grid, x, -1 + 2 * x)

df_y = data.frame(x = y_grid, y = f_yx)
ggplot(df_y, aes(x, y)) + geom_point(size = 0.9) +
    scale_x_continuous(breaks = seq(-1.5, 1.5, by = 0.5))

df_yx = data.frame(x = X0[,2], y = y0)
ggplot(df_yx, aes(x, y)) + geom_point()



# second simulation ------------------------------------------------------------