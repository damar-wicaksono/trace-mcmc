#' Calculate the variance matrix due to kriging prediction in the LMC model
#'
#' \code{CalcKrigingVarMat} basically reconstructs the variance matrix in the  
#' original output space due to the standard deviation of the prediction 
#' made on the Principal component scores by GP metamodel.
#' LMC is a linear model with coefficients of the model is predicted by GP 
#' model. 
#' So you can think of it as uncertainty propagation.
#'
#' @param pc_scores_sd Numeric vector. The standard deviation of the kriging
#'  prediction of the principal component scores
#' @param trc_pca_lds Matrix. The principal component loadings, each column is
#'  a principal component loading
#' @param idx Numeric vector. The select indices of the full output map
#' @return Matrix. The variance at the output space at select points due to 
#'  kriging variance
CalcKrigingVarMat <- function(pc_scores_sd, trc_pca_lds, idx)
{
    num_pc <- length(pc_scores_sd)
    pc_scores_var <- diag(pc_scores_sd^2)

    kriging_var_mat <- trc_pca_lds[idx, 1:num_pc] %*% pc_scores_var %*%
        t(trc_pca_lds[idx, 1:num_pc])
    
    return(kriging_var_mat)
}
