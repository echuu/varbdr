


# (D x max_iter) matrix for each cluster; contains each iter's value for
# the mean vector
# m_k_hist  = array(0, c(D, max_iter, K))
# mu_k_hist = array(0, c(D, max_iter, K)) 

# store frob norm of matrix per iteration for each cluster
# V_k_hist = matrix(0, K, max_iter)
# Q_k_hist = matrix(0, K, max_iter)  


# plot differences in Q_k, V_k

V_k_conv = t(theta1_2$V_k_hist[,1:theta1_2$curr-1]) # curr x K
Q_k_conv = t(theta1_2$Q_k_hist[,1:theta1_2$curr-1]) # curr x K

V_k_conv_long = melt(V_k_conv)
Q_k_conv_long = melt(Q_k_conv)

names(V_k_conv_long) = c("iter", "cluster", "value")
names(Q_k_conv_long) = c("iter", "cluster", "value")


ggplot(V_k_conv_long, aes(x = iter, y = value, col = factor(cluster))) + 
    geom_point(size = 1) + theme_bw() + ggtitle("frob norm for prec. matrix Vk")

ggplot(Q_k_conv_long, aes(x = iter, y = value, col = factor(cluster))) + 
    geom_point(size = 1) + theme_bw() + ggtitle("frob norm for prec. matrix Qk")


# plot differences in m_k, mu_k
m_1_conv = theta1_2$m_k_hist[,,1][1:theta1_2$curr-1]
m_2_conv = theta1_2$m_k_hist[,,2][1:theta1_2$curr-1]

mu_1_conv = theta1_2$mu_k_hist[,,1][1:theta1_2$curr-1]
mu_2_conv = theta1_2$mu_k_hist[,,2][1:theta1_2$curr-1]


m_df = data.frame(m_1 = m_1_conv, m_2 = m_2_conv)
m_long = melt(m_df)
m_long$iter = rep(1:nrow(m_df), 2)
names(m_long) = c("cluster", "value", "iter")


mu_df = data.frame(mu_1 = mu_1_conv, mu_2 = mu_2_conv)
mu_long = melt(mu_df)
mu_long$iter = rep(1:nrow(mu_df), 2)
names(mu_long) = c("cluster", "value", "iter")

ggplot(m_long, aes(x = iter, y = value, col = factor(cluster))) + 
    geom_point(size = 1) + theme_bw() + ggtitle("values of each component m_k")

ggplot(mu_long, aes(x = iter, y = value, col = factor(cluster))) + 
    geom_point(size = 1) + theme_bw() + ggtitle("values of each component mu_k")



