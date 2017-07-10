#' Helper function to get the vector of indices to subset the full TC output
#'
#' Stationary model of bias only applies to the region of the transient before
#' quenching, \code{GetTimeIdxTC} get the indices of the full model output 
#' prior to quenching in all axial locations.
#'
#' @param trc_data_bias List. List of compiled results from bias model 
#'  training runs
#' @return Numeric vector. The indices from the full model output before 
#'  quenching where stationary model applies, to be compared with experimental
#'  data
GetTimeIdxTC <- function(trc_data_bias)
{
    exp_data <- trc_data_bias$exp_data[[1]]
    trc_time <- trc_data_bias$time
    trc_nom  <- trc_data_bias$nominal[,1:8]
    trc_runs <- trc_data_bias$replicates[,,1:8]
    
    res_time <- GetTimeExpTC(exp_data, trc_time, trc_nom, trc_runs)
    time_idx <- FlattenTimeIdxExp(
        res_time$trc_idx, length(trc_time), c(1)
    )
    
    return(time_idx)
}