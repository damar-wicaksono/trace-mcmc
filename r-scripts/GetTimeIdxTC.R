#' Helper function to get the vector of indices to subset the full TC output
#'
#' Stationary model of bias only applies to the region of the transient before
#' quenching, \code{GetTimeIdxTC} get the indices of the full model output 
#' prior to quenching in all axial locations.
#'
#' @param trc_data_bias List. List of compiled results from bias GP model 
#'  training runs
#' @param excl_idx Numeric vector. Indices to be excluded 
#'  (default = initial time point for each axial location)
#' @param d_t Numeric. The time step size in the uniform time grid of full 
#'  output space (default = 0.1 [s])
#' @param temp_min Numeric. The minimum temperature in [K] where the rod can  
#'  be assumed to be quenched (default = 400.0 [K])
#' @return Numeric vector. The indices from the full model output before 
#'  quenching where stationary model applies, to be compared with experimental
#'  data
GetTimeIdxTC <- function(trc_data_bias, 
                         excl_idx = c(1), d_t = 0.1, temp_min = 400.)
{
    exp_data <- trc_data_bias$exp_data[[1]]
    trc_time <- trc_data_bias$time
    trc_nom  <- trc_data_bias$nominal[,1:8]
    trc_runs <- trc_data_bias$replicates[,,1:8]
    
    res_time <- GetTimeExpTC(
        exp_data,
        trc_time,
        trc_nom,
        trc_runs,
        dt,
        temp_min)

    time_idx <- FlattenTimeIdxExp(
        res_time$trc_idx,
        length(trc_time),
        excl_idx
    )
    
    return(time_idx)
}
