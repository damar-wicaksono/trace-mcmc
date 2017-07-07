#' Calculate the variance matrix due to bias model in the LMC model
#'
#' \code{CalcTruncationVarMat} reconstructs the variance matrix in the original 
#' output space due to the bias modelled as a Gaussian process.
#'
#' @param xx_df Dataframe.
#' @param xx_var Numeric.
#' @param trc_dis_cov CovTensorProduct (DiceKriging).
#' @return Matrix. Variance matrix due to the bias modeled as a Gaussian 
#'  process evaluated at the full output space
#'
CalcBiasVarMat <- function(xx_df, xx_var, trc_dis_cov)
{
    cov_bias <- trc_dis_cov
    # Override the variance parameter
    cov_bias@sd2 <- xx_var
    # Compute the bias variance
    bias_var_mat <- covMatrix(cov_bias, as.matrix(xx_df))$C

    return(bias_var_mat)
}
