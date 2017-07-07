#' Compute the Log-Likelihood of Normal Model
#'
#' The inversion of the variance matrix is done by SVD
#'
#' @param exp_vec Numeric vector. Vector of experimental data
#' @param ave_vec Numeric vector. Vector of the prediction mean
#' @param var_mat Matrix. Variance matrix
#' @return Numeric. The log-likelihood
CalcLogLikelihood <- function(exp_vec, ave_vec, var_mat)
{
    # Variance Matrix Inversion
    svd_var_mat <- svd(var_mat)
    inv_var_mat <- svd_var_mat$v %*% diag(1/svd_var_mat$d) %*% t(svd_var_mat$u)
    
    k <- length(exp_vec)
    b1 <- matrix(exp_vec - ave_vec, ncol = k) %*% inv_var_mat %*% 
        matrix(exp_vec - ave_vec, nrow = k)
    b2 <- determinant(var_mat, logarithm = T)$modulus[1]
    b3 <- k * log(2 * pi)
    
    return(-0.5 * (b1 + b2 + b3))
}
