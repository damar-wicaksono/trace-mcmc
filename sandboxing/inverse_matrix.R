full_output_dim <- dim(var_mat_disc_tc)[1] + dim(var_mat_discrep_dp)[1] + dim(var_mat_discrep_co)[1]
ave_vec_full <- c(ave_vec_tc, ave_vec_dp, ave_vec_co)
var_mat_full <- matrix(0, nrow = full_output_dim, ncol = full_output_dim)

var_mat_full[1:length(ave_vec_tc),1:length(ave_vec_tc)] <- var_mat_krig_tc + var_mat_trun_tc + var_mat_disc_tc
var_mat_full[(length(ave_vec_tc)+1):(length(ave_vec_tc)+length(ave_vec_dp)),(length(ave_vec_tc)+1):(length(ave_vec_tc)+length(ave_vec_dp))] <- var_mat_kriging_dp + var_mat_truncat_dp + var_mat_discrep_dp
var_mat_full[(length(ave_vec_tc)+length(ave_vec_dp)+1):(length(ave_vec_tc)+length(ave_vec_dp)+length(ave_vec_co)),(length(ave_vec_tc)+length(ave_vec_dp)+1):(length(ave_vec_tc)+length(ave_vec_dp)+length(ave_vec_co))] <- var_mat_kriging_co + var_mat_truncat_co + var_mat_discrep_co

exp_vec_full <- c(tc_exp_vector, dp_exp_vector, co_exp_vector)
log(
    dmvnorm(exp_vec_full, 
            mean = ave_vec_full, 
            sigma = var_mat_full)
)

log(
    dmvnorm(exp_vec_full[128:199], 
            mean = ave_vec_full[128:199], 
            sigma = var_mat_full[128:199,128:199])
)
log(
    dmvnorm(exp_vec_full[1:127], 
            mean = ave_vec_full[1:127], 
            sigma = var_mat_full[1:127,1:127])
)
log(
    dmvnorm(exp_vec_full[200:206], 
            mean = ave_vec_full[200:206], 
            sigma = var_mat_full[200:206,200:206])
)
-572.9385 -583.7696-8.928393

length(ave_vec_full)
full_output_dim
exp_vec_full
127 + 72

solve(var_mat_full[1:127,1:127])
log(
    dmvnorm(exp_vec_full[1:199], 
            mean = ave_vec_full[1:199], 
            sigma = var_mat_full[1:199,1:199])
)

1 / det(var_mat_full[1:127,1:127])^0.5

L <- chol(var_mat_full)
Q <- t(L)
determinant(L, logarithm = T)$modulus[1]

diag(t(solve(L)) %*% solve(L))
diag(solve(var_mat_full))

svd_var_mat_krig_tc <- svd(var_mat_trun_tc)
diag(var_mat_trun_tc %*% svd_var_mat_krig_tc$v %*% 
         diag(1/svd_var_mat_krig_tc$d) %*% t(svd_var_mat_krig_tc$u))

det(Q)
determinant(var_mat_full, logarithm = T)

a1 <- -0.5 * matrix(tc_exp_vector - ave_vec_tc, ncol = 127) %*% solve(var_mat_full[1:127,1:127]) %*% matrix(tc_exp_vector - ave_vec_tc, nrow = 127)
a2 <- -0.5 * determinant(var_mat_full[1:127,1:127], logarithm = T)$modulus[1]
a3 <- -0.5 * 127 * log(2*pi)

a1 + a2 + a3


b1 <- -0.5 * matrix(exp_vec_full - ave_vec_full, ncol = 206) %*% chol2inv(chol(var_mat_full)) %*% matrix(exp_vec_full - ave_vec_full, nrow = 206)
b1 <- -0.5 * matrix(exp_vec_full - ave_vec_full, ncol = 206) %*% solve(var_mat_full) %*% matrix(exp_vec_full - ave_vec_full, nrow = 206)
diag(solve(var_mat_full))
diag(chol2inv(chol(var_mat_full)))
b2 <- -0.5 * determinant(var_mat_full, logarithm = T)$modulus[1]
b3 <- -0.5 * 206 * log(2*pi)
b1 + b2 + b3
b1

L <- chol(var_mat_full[1:127,1:127])
Q <- t(L)
determinant(L, logarithm = T)$modulus[1]

diag(chol2inv(L))
diag(solve(var_mat_full[1:127,1:127]))

svd_var_mat_full <- svd(var_mat_full[1:127,1:127])
diag(svd_var_mat_full$v %*% diag(1/svd_var_mat_full$d) %*% t(svd_var_mat_full$u))


# Full Covariance Matrix
diag(chol2inv(chol(var_mat_full)))
diag(solve(var_mat_full))
svd_var_mat_full <- svd(var_mat_full)
diag(svd_var_mat_full$v %*% diag(1/svd_var_mat_full$d) %*% t(svd_var_mat_full$u))

diag(var_mat_full %*% chol2inv(chol(var_mat_full)))
diag(var_mat_full %*% solve(var_mat_full))
diag(var_mat_full %*% svd_var_mat_full$v %*% diag(1/svd_var_mat_full$d) %*% t(svd_var_mat_full$u))


# Computing Inverse Matrix using R, battery included
# There are, say, 4 ways of computing inverse:
mat_test <- matrix(c(5,1,1,3),2,2)
# 1. Direct
mat_test %*% solve(mat_test) # Look's good enough
# 2. Cholesky, indirect
L <- chol(mat_test)
mat_test %*% t(solve(L)) %*% solve(L) # Not so good
# 3. Cholesky, direct
mat_test %*% chol2inv(L) # Pretty good 
# 4. SVD
svd_mat_test <- svd(mat_test)
mat_test %*% svd_mat_test$v %*% diag(1/svd_mat_test$d) %*% t(svd_mat_test$u) # Again pretty good

