#' Get restricted time index for temperature output (TC <= quenching)
#'
#' The bias modeling used in this study is only valid for part of the transient 
#' that can be modelled using stationary gaussian process.
#' The quenching of the rod wall is thus problematic.
#' A work around of this problem is simply to compare the prediction and 
#' experimental data where quenching has not occured.
#' This approach can be justified due to the fact that increasing the agreement 
#' between the prediction and measurement data before the quenching is enough
#' to constraint the most influential parameters.
#' Indeed, parameters whose impact is more localized in time risks not to be
#  well-calibrated. 
#' 
#' @param exp_data Matrix. Time-Experimental Data pairs, 1st column is the time
#'  column 2-9 are the temperature measurement data in [degC] at different 
#'  axial locations (2 is the bottom while 9 is the top).
#' @param trc_time Numeric vector. TRACE output time points
#' @param trc_nom Matrix. TRACE temperature output, nominal run
#' @param trc_runs Array. TRACE temperature output, bias model training runs
#' @param d_t Numeric. TRACE time step size
#' @param temp_min Numeric. Minimum temperature to be considered "quenched"
#' @return List. $exp_idx is the number of experimental points before 
#'  quenching, $trc_idx is the list with time indices for each axial 
#'  location, while $trc_pts is the list with time points for each axial 
#'  location, taken from the experimental data table with temperature 
#'  measurement up to before quenching occur for all the training runs 
#'  for the bias model
GetTimeExpTC <- function(exp_data, trc_time, trc_nom, trc_runs, d_t, temp_min)
{
    n_samples <- dim(trc_runs)[2]  # Number of training samples
    exp_time  <- exp_data[, 1]     # Experimental time points
    exp_data  <- exp_data[,-1]     # Experimental TC1-TC8

    # Pre-process the Data --------------------------------------------------------
    # Read the experimental time of quenching, 
    # based on maximum difference between two consecutive points
    diff_exp_tc <- exp_data[1:(nrow(exp_data)-1),1:8] - 
        exp_data[2:nrow(exp_data),1:8]
    exp_tquench_idx <- apply(diff_exp_tc, 2, which.max)

    # Loop over 8 axial locations
    trc_time_restricted_pts <- list()
    trc_time_restricted_idx <- list()
    for (j in 1:dim(trc_nom)[2])
    {
        # Compute the time of quenching at an axial location for all samples
        trc_runs_tquench <- c()
        for (i in 1:n_samples)
        {
            trc_runs_tquench <- c(trc_runs_tquench,
                                  kneedle(trc_time, trc_runs[,i,j]))
        }
    
        # Compute the nominal time of quenching
        trc_nom_tquench <- c(trc_runs_tquench, kneedle(trc_time, trc_nom[,j]))
    
        # The index of time step for TRACE output up to quenching
        trc_time_idx <- exp_time[1:exp_tquench_idx[j]] / d_t + 1 
    
        # Loop over backward from the end of transient
        k <- 0
        trc_exp_time <- trc_time[trc_time_idx]
        for (i in length(trc_time_idx):1)
        {
            # Remove one point of experimental data if:
            if (min(trc_runs_tquench) < trc_exp_time[i])  
            {
                # 1. The minimum time of quenching across bias training runs is 
                #    smaller than a given location of time of of quenching 
                #    experimental time
                k <- k + 1
            } else if (min(trc_runs[trc_time_idx[i],,j]) %/% 100. <= temp_min %/% 100.)
            {
                # 2. The minimum temperature of any replication is less 
                #    or equal temp_min [K]
                k <- k + 1
            } else
            {
                # Otherwise exit the loop, go to next axial location
                break
            }
        }
    
        # Exclude Data Points where early quenching happens
        exp_tquench_idx[j] <- exp_tquench_idx[j] - k
        trc_time_idx <- trc_time_idx[1:i]
        trc_time_restricted_idx[[j]] <- trc_time_idx
        trc_time_restricted_pts[[j]] <- trc_time[trc_time_idx]
    }

    output <- list(
        exp_idx = exp_tquench_idx,
        trc_idx = trc_time_restricted_idx,
        trc_pts = trc_time_restricted_pts
    )

    return(output)
}
