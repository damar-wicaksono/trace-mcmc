#' Helper function to get the vector of TC experimental data points
#'
#' The measurement data is directly converted to [K] from [degC].
#'
#' @param trc_data_bias List. List of compiled results from TRACE bias model 
#' training runs
#' @param d_t Numeric. The time step size in the uniform time grid of full 
#'  output space (default = 0.1 [s])
#' @param temp_min Numeric. The minimum temperature in [K] where the rod can  
#'  be assumed to be quenched (default = 400.0 [K])
#' @return Numeric vector. The TC experimental data points selected before 
#'  quenching
GetExpDataTC <- function(trc_data_bias, d_t = 0.1, temp_min = 400.)
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
    
    exp_vec <- c()
    for (j in 1:8)
    {
        exp_vec <- c(exp_vec, exp_data[2:res_time$exp_idx[j], j+1])
    }
    exp_vec <- exp_vec + 273.15 # Convert to [Kelvin]
    
    return(exp_vec)
}
