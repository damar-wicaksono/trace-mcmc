#' Calculate the average/mean vector in the LMC model
#' 
#' \code{CalcAveVec} basically reconstructs the original output space based 
#' on the prediction made on the Principal component scores by GP metamodel.
#'
#' @param pc_scores_mu Numeric vector. The mean of GP prediction of PC scores
#' @param trc_pca_ave Numeric vector. Mean of full output space from training
#' @param trc_pca_lds Matrix. PC Loadings from PCA training
#' @param idx Numeric vector. The select indices of the full output map
#' @return Numeric vector. The reconstructed outputs based on predicted PC 
#'  scores and principal componets, at select indices
CalcAveVec <- function(pc_scores_mu, trc_pca_ave, trc_pca_lds, idx)
{
    # Matrix operation to reconstruct the full output space
    pc_scores_mu <- matrix(pc_scores_mu, nrow = length(pc_scores_mu))
    
    ave_vec <- trc_pca_ave + trc_pca_lds[, 1:length(pc_scores_mu)] %*% pred_mu
    
    # Return only select points in the full output space
    ave_vec <- ave_vec[idx]
    
    return(ave_vec)
}
