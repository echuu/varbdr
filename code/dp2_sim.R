

setwd(getwd())

source("globals.R")

setwd(HOME_DIR)
source(DP_BDR) 
source(DENSITY)


N = 500
K = 2
synth_data_1d = r_dpmix2(N)

y      = synth_data_1d$y
X      = synth_data_1d$X
params = synth_data_1d$param_mat
y_grid = seq(-0.5, 1.5, length.out = 500)

# 6th sublot in Figure 3
df_yx = data.frame(x = X[,2], y = y)
ggplot(df_yx, aes(x = x, y = y)) + geom_point() + stat_smooth() + theme_bw()


# 1st subplot

x         = c(0.15, 0.25, 0.49, 0.75, 0.88, 0.95)
f_yx      = matrix(0, nrow = N, ncol = length(x))

for(i in 1:length(x)) {
    # f_yx[,i] = d_dpmix2(y_grid, x[i], -1 + 2 * x[i])
    f_yx[,i] = d_dpmix2(y_grid = y_grid, x = c(1, x[i]), 
                        params = c(x[i], x[i]^4))
}

df_y      = data.frame(x = y_grid, f_yx)
df_y_long = melt(df_y, id.vars = "x")
levels(df_y_long$variable) = c("x = 0.15", "x = 0.25", "x = 0.49",
                               "x = 0.75", "x = 0.88", "x = 0.95")
ggplot(df_y_long, aes(x, y = value)) + geom_point(size = 0.9) +
    scale_x_continuous(breaks = seq(-0.5, 1.5, by = 0.5)) + 
    facet_wrap(~variable)


df_yx = data.frame(x = X[,2], y = y)
ggplot(df_yx, aes(x, y)) + geom_point() + ggtitle("response vs. covariate") + 
    stat_smooth(se = FALSE) + theme_bw()


# run cavi ---------------------------------------------------------------------

# (1.1) covariate-INDEPENENT vb
source(paste(COV_INDEP, VARBDR, sep = '/'))      # load cov-indep varbdr.R
theta1_1 = varbdr(y = y, X = X, K)               # store CAVI results

# (1.2) covariate-DEPENDENT vb
source(paste(COV_DEP, VARBDR, sep = '/'))        # load cov-dep varbdr.R
theta1_2 = varbdr(y = y, X = X, K)               # store CAVI results


# save variational parameters
# saveRDS(theta1_1, file = "dp_results/theta_indep_dp2_K2.RDS")
# saveRDS(theta1_2, file = "dp_results/theta_dep_dp2_K2.RDS")

plotELBO(theta1_1)
plotELBO(theta1_2)


# generate plots ---------------------------------------------------------------

theta1      = list(theta1_1, theta1_2)      # list of var. params for each alg
DATA_GEN_ID = DP_MIX2

# generate plot overlays of the percentiles
source(DENSITY)
x = c(0.10, 0.25, 0.49, 0.75, 0.88, 0.95)
approx_d  = list(py_0, py_bouch)            # list of approx density functions
den_label = c("cov-indep", "cov-dep (b)")   # labels for each approx density
p_list1 = xQuantileDensity(x, params, d_dpmix2,
                           approx_d, den_label, theta1, K,
                           DP_MIX2,
                           y_grid,
                           DP2 = TRUE)

# p_list1[[3]]
multiplot(plotlist = p_list1, cols = 3)





