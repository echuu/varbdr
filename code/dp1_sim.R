
# dp1_sim.R() -- conditional density estimation, true density coming from 
#                the first numerical example in Dunson & Park paper

setwd(getwd())

source("globals.R")

setwd(HOME_DIR)
source(DP_BDR) 
source(DENSITY)

## generate data from true density ---------------------------------------------
N = 500
K = 4
dp1_synth = r_dpmix1(N)
y0        = dp1_synth$y
X0        = dp1_synth$X
mu1       = dp1_synth$param_vec
y_grid    = seq(-1.5, 1.5, length.out = 1000)

## generate figures for conditional density for increasing values of x -----
x         = c(0.15, 0.25, 0.49, 0.75, 0.88, 0.95)
f_yx      = matrix(0, nrow = length(y_grid), ncol = length(x))
for(i in 1:length(x)) {
    f_yx[,i]      = d_dpmix1(y_grid, x[i], -1 + 2 * x[i])
}

df_y      = data.frame(x = y_grid, f_yx)

df_y_long = melt(df_y, id.vars = "x")
levels(df_y_long$variable) = c("x = 0.15", "x = 0.25", "x = 0.49",
                               "x = 0.75", "x = 0.88", "x = 0.95")
ggplot(df_y_long, aes(x, y = value)) + geom_point(size = 0.9) +
    scale_x_continuous(breaks = seq(min(y_grid), max(y_grid), by = 1)) + 
    facet_wrap(~variable)

df_yx = data.frame(x = X0[,2], y = y0)
ggplot(df_yx, aes(x, y)) + geom_point() + ggtitle("response vs. covariate")


## run vb algorithm to find optimal parameters ---------------------------------

# (1.1) covariate-INDEPENENT vb
source(paste(COV_INDEP, VARBDR, sep = '/'))      # load cov-indep varbdr.R
theta1_1 = varbdr(y = y0, X = X0, K)             # store CAVI results

# (1.2) covariate-DEPENDENT vb
source(paste(COV_DEP, VARBDR, sep = '/'))        # load cov-dep varbdr.R
theta1_2 = varbdr(y = y0, X = X0, K)             # store CAVI results


# save variational parameters
# saveRDS(theta1_1, file = "dp_results/theta_indep_dp1_K4.RDS")
# saveRDS(theta1_2, file = "dp_results/theta_dep_dp1_K4.RDS")

theta1_1   = readRDS(file = "dp_results/theta_indep_dp1_K4.RDS")
old_theta2 = readRDS(file = "dp_results/theta_dep_dp1_K4.RDS")



# plot ELBO
plotELBO(theta1_1)
plotELBO(theta1_2)


# generate plots ---------------------------------------------------------------
source(DENSITY)
source(paste(COV_DEP, VARBDR, sep = '/'))
theta1      = list(theta1_2)      # list of var. params for each alg
DATA_GEN_ID = DP_MIX1

# generate plot overlays of the percentiles
x = c(0.15, 0.25, 0.49, 0.75, 0.88, 0.95)
approx_d  = list(py_bouch)        # list of approx density functions
den_label = c("cov-dep")          # labels for each approx density
params = c(-1 + 2 * x)
p_list1 = xQuantileDensity(x, params, d_dpmix1,
                           approx_d, den_label, theta1, K,
                           DATA_GEN_ID,
                           y_grid)
multiplot(plotlist = as.list(sapply(p_list1, function(x) x[1])), cols = 3)

# look at density curve over iterations for cov-dependent density --------------

n = 100 # row that has 2nd element ~ 0.15 -> compare to p_list1[[1]] above
mu1 = -1 + 2 * X0[n,2] # mean for y | x_n

plotELBO(theta1_2)

iters = c(2, 3, 4, 5, 8, 20, 23, 24, 25, 26, 27, 28, 40, 400)
plot_df = matrix(0, nrow = N, ncol = length(iters))
for (i in 1:length(iters)) {
    plot_df[,i] = theta1_2$dc[[iters[i]]][n,]
}
plot_df = data.frame(plot_df)
names(plot_df) = paste("n", iters, sep = '_')

plot_df$y = d_dpmix1(y_grid, X0[n,], mu1)
plot_df$x = y_grid
df_long = melt(plot_df, id.vars = "x")

ggplot(df_long, aes(x = x, y = value)) + 
    geom_line(aes(x = df_long$x, y = df_long$value, 
                  col = df_long$variable), size = 1.0) + theme_bw()

# by 40 iterations, the approx is already same as iteration 400


# compare with np results ------------------------------------------------------
library("np")

X = data.frame(X0[,2])
y = data.frame(y0)
xy_df = data.frame(y = y, x = X)
names(xy_df) = c("x", "y")

fy.x = npcdens(y ~ x, xy_df)

overlays = vector("list", length(x))

for (i in 1:length(x)) {
    y_eval = data.frame(y = y_grid, x = x[i])
    fy_eval = predict(fy.x, newdata = y_eval)
    
    fy_df = data.frame(y = y_grid, np = fy_eval)
    
    # ggplot(fy_df, aes(x = y, y = np)) + geom_line()
    
    # overlay np estimate with the other approximations ------------------------
    
    # x = 0.15
    full_approx_df = merge(p_list1[[i]]$approx_df, fy_df, by = "y")
    approx_df_long = melt(full_approx_df, 
                          measure.vars = c(den_label, "true", "np"))
    
    overlays[[i]] = ggplot(approx_df_long, 
                           aes(x = y, y = value, colour = variable)) + 
        geom_line(size = 0.8) + labs(x = "y", y = "p(y)") + theme_bw() +
        theme(legend.position = "none")
}

multiplot(plotlist = overlays, cols = 3)

# dp1_sim.R
