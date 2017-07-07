#' Calculate the mean and standard deviation of PC scores from GP metamodel
#'
#' The implementation avoid using \code{for} to loop over the specified number
#' of \code{km} objects.
#' This is achieved by using \code{mapply()}.
#' \code{mapply()} will produce 2-dimensional matrix of lists.
#' The row is the list of an output from calling predict for all km objects
#' The column is the list of all predict outputs for a given km object
#'
#' @param xx_df Dataframe. The input parameter values as a dataframe, requires
#'  a consistent number of columns and naming as the ones in the km objects.
#' @param trc_gps List of km objects (DiceKriging). The trained kriging model 
#'  for each of the PC scores.
#' @param num_pc Numeric. The number of PC used in the reconstruction.
CalcPCScores <- function(xx_df, trc_gps, num_pc)
{
    pc_scores <- mapply(predict, 
                        object = trc_gps[c(1:num_pc)], 
                        newdata = replicate(num_pc, xx_df, simplify = FALSE),
                        type = "SK")

    return(list(mu = unlist(pc_scores["mean",]),
                sd = unlist(pc_scores["sd",])))
}
