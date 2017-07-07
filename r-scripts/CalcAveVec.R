#' Calculate the average/mean vector of the LMC model
#' 
#' \code{CalcAveVec} basically reconstructs the original output space based 
#' on the prediction made on the Principal component scores.
#' 
#' The implementation avoid using \code{for} to loop over the specified number
#' of \code{km} objects.
#' This is achieved by using \code{mapply()}.
#' \code{mapply()} will produce 2-dimensional matrix of lists.
#' The row is the list of an output from calling predict for all km objects
#' The column is the list of all predict outputs for a given km object
#'
#' @param xx_df Dataframe. The input parameter values as a dataframe, requires
#'  a consistent number of columns and naming as the ones in the km objects
#' @param trc_gps List of km objects (DiceKriging). The trained kriging model 
#'  for each of the PC scores
#' @param trc_pca_ave Numeric vector. Mean of full output space from training
#' @param trc_pca_lds Matrix. PC Loadings from PCA training
#' @param num_pc Numeric. The number of PC used in the reconstruction
#' @param idx Numeric vector. The select indices of the full output map.
#' @return Numeric vector. The reconstructed outputs based on predicted PC 
#'  scores and principal componets, at select indices.
#'
CalcAveVec <- function(xx_df, trc_gps, 
                       trc_pca_ave, trc_pca_lds, num_pc, idx)
{
    pred_mu <- c()
    pred_mu <- mapply(predict, 
                      object = trc_gps[c(1:num_pc)], 
                      newdata = replicate(num_pc, xx_df, simplify = FALSE),
                      type = "SK")

    # Flatten the list of "mean" output from predict()
    pred_mu <- unlist(pred_mu["mean",])
    
    # Matrix operation to reconstruct the full output space
    pred_mu <- matrix(pred_mu, nrow = num_pc)
    
    ave_vec <- trc_pca_ave + trc_pca_lds[, 1:num_pc] %*% pred_mu
    
    # Return only select points in the full output space
    ave_vec <- ave_vec[idx]
    
    return(ave_vec)
}
