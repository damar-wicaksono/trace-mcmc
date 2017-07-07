#' Calculate the variance matrix due to kriging prediction in the LMC model
#'
#' \code{CalcTruncationVarMat} basically reconstructs the variance matrix in 
#  the original output space due to the truncation of the principal components.
#'
#' The truncation error coefficient is modelled as a normal distribution with 
#' variance 1 (because the PC loading has the magnitude already).
#' As such the variance on the output space due to the truncation is simply a 
#' square of the principal component loadings not used in the reconstruction of
#' the output space.
#'
#' @param trc_pca_lds Matrix. The principal component loadings, each column is
#'  a principal component loading
#' @param num_pc Numeric. The number of PC components used to reconstruct the 
#'  output space
#' @param idx Numeric vector. The select indices of the full output map
#' @return Matrix. The variance at the output space at select points due to
#'  truncation error
CalcTruncationVarMat <- function(trc_pca_lds, num_pc, idx)
{
    num_pc_total <- dim(trc_pca_lds)[2]

    # PC Truncation Variance
    pca_lds <- trc_pca_lds[idx, (num_pc + 1):num_pc_total]
    truncation_var_mat <- pca_lds %*%  
        diag(replicate(1, num_pc_total - num_pc)) %*% t(pca_lds)
    
    return(truncation_var_mat)
}
