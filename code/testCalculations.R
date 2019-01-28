



# re-visit m-step calculation

N = 10
K = 4
r_nk = runif(N)
r_diag = diag(r_nk) / sum(r_nk)


X = matrix(rnorm(N * K), nrow = 10)

vec_prod = t(X) %*% sqrt(r_diag) %*% X

mat_sum = 0
for (n in 1:N) {
    mat_sum = mat_sum + r_nk[n] * crossprod(t(X[n,]))
}
mat_sum

mat_sum = 0
for (n in 1:N) {
    mat_sum = mat_sum + crossprod(t(X[n,]))
}

y = rnorm(N)
ry_x = 0
for (n in 1:N) {
    ry_x = ry_x + r_nk[n] * y[n] * X[n,]
}

t(X) %*% (r_nk * y)
