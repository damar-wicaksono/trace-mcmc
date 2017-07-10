#' Helper function to get the total number of TC exp. points before quenching
#'
#' The measurement data is directly converted to [K] from [degC].
#'
#' @param trc_data_bias List. List of compiled results from TRACE bias model 
#' training runs
#' @param d_t Numeric. The time step size in the uniform time grid of full 
#'  output space (default = 0.1 [s])
#' @param temp_min Numeric. The minimum temperature in [K] where the rod can  
#'  be assumed to be quenched (default = 400.0 [K])
#' @return Numeric vector. The number of TC experimental data points per axial
#'  location before quenching
GetExpTimeQuenchIdx <- function(trc_data_bias, d_t = 0.1, temp_min = 400.)
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
        d_t,
        temp_min)

    return(res_time$exp_idx)
}
