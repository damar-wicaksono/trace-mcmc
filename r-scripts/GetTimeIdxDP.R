#' Helper function to get the vector of indices to subset the full DP output
#'
#' The function is straightforward as all points available in the experimental
#' dataset can be compared with the prediction.
#' There were 4 pressure drop measurements and all had the same time points.
#'
#' @param trc_data_bias List. List of compiled results from bias model 
#'  training runs
#' @param d_t Numeric. The time step size in the uniform time grid of full 
#'  output space (default = 0.1 [s])
#' @return Numeric vector. The indices from the full model output for all 
#'  pressure drop measurements.
GetTimeIdxDP <- function(trc_data_bias, d_t = 0.1)
{
    exp_data <- trc_data_bias$exp_data[[2]]
    res_time <- GetTimeIdxTrc(exp_data[, 1], d_t = d_t)
    res_time <- list(res_time, res_time, res_time, res_time)
    time_idx <- FlattenTimeIdxExp(res_time, length(trc_time))

    return(time_idx)
}
