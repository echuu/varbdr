

N = 10
K = 3
D = 3

# variables needed: r_nk, tau_k, X, beta
# X (N x K)
set.seed(123)
X = matrix(rnorm(N * K), N, K)

# beta (D x K)
beta = matrix(rnorm(D * K), D, K)

# r_nk (N x K)
r_nk = matrix(rnorm(N * K), N, K)

# tau_k (K x 1)
tau_k = rnorm(K)

sum_res = 0

for (d in 1:D) {
    for (k in 1:K) {
        for (n in 1:N) {
            for (j in 1:D) {
                
                if (j == d) {
                    next
                }
                
                sum_res = sum_res + r_nk[n,k] * tau_k[k] * X[n,d] * X[n,j] * 
                    beta[d,k] * beta[j,k]
                
            }
        }
    }
}

sum_res # -2.56325

sum_res = 0

for (d in 1:D) {
    
    j_sum = 0
    
    for (j in 1:D) {
        
        if (j == d) {
            next
        }
        
        R_dj = matrix(0, K, K)
        
        # construct the (K x K) R_dj matrix
        for (k in 1:K) {
            R_dj[k,k] = tau_k[k] * sum(r_nk[,k] * X[,j] * X[,d])
        }     
        
        # compute the summation over j != d
        j_sum = j_sum + R_dj %*% beta[j,]

    }
    
    # compute the summation over d
    sum_res = sum_res + t(beta[d,]) %*% j_sum 
}

sum_res # -2.56325



d = 1
j = 2

t(beta[d,]) %*% diag(tau_k * r_nk[1,]) %*% beta[j,]

test_sum = 0
for (k in 1:K) {
    test_sum = test_sum + tau_k[k] * r_nk[1,k] * beta[d,k] * beta[j,k]
}
test_sum



vec_res = 0
for (d in 1:D) {
    
    for (j in 1:D) {
        if (j == d) {
            next
        }
        
        
        j_sum = 0
        
        
        for (n in 1:N) {
            R_n = diag(tau_k * r_nk[n,])
            
            n_mat = X[n,d] * X[n,j] * R_n
        }
        
        j_sum = n_mat %*% beta[j,]
        
        
    }
    
    vec_res = vec_res + t(beta[d,]) %*% j_sum
    
}
vec_res


vec_res = 0
for (d in 1:D) {
    
    j_sum = 0
    for (j in 1:D) {
        if (j == d) {
            next
        }
        
        zeta_dj = diag(0, K)
        x_prod = X[,d] * X[,j]
        for (k in 1:K) {
            zeta_dj[k,k] = tau_k[k] * sum(r_nk[,k] * x_prod)
        }
        
        j_sum = j_sum + zeta_dj %*% beta[j,]
    }
    
    vec_res = vec_res + t(beta[d,]) %*% j_sum
    
    
}

vec_res












