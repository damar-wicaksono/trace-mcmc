#' Set initial values of the chain in the walker
#'
#' @param num_model_params Numeric. The number of model parameters
#' @param init_sd_bias List. The iniitial bias model scale parameter(s)
#' @param n_walks Numeric. Number of walkers
#' @param n_iters Numeric. Number of iterations
#' @param trc_outputs Character vector. The type of outputs to be considered
#' @return 3D array, 1st dimension is the container for each parameter, 
#'  2nd dimension is the container for each walker,
#'  3rd dimension is the sample per each iteration
setInitialValues <- function(num_model_params, 
                             init_sd_bias, n_walks, n_iters, trc_outputs)
{
    sd_multiplier <- list("tc" = 10., "dp" = 100., "co" = 0.1)
    num_params <- num_model_params + length(trc_outputs)

    # Set up initial values
    if (length(trc_outputs) < 3)
    {
        # All outputs
        post_samples <- array(NA, dim = c(num_params, n_walks, n_iters))
        for (i in 1:num_model_params)
        {
            # Initial values at nominal then perturbed
            post_samples[i,,1] = 0.5 + 1e-4 * rnorm(n_walks)
        }
        # The last parameter is the scale parameter of the bias model
        post_samples[num_params,,1] = init_sd_bias[[trc_outputs[1]]] + 
            sd_multiplier[[trc_outputs[1]]] * rnorm(n_walks)
    } else
    {
        post_samples <- array(NA, dim = c(num_params, n_walks, n_iters))
        for (i in 1:num_model_params)
        {
            # Initial values at nominal then perturbed
            post_samples[i,,1] = 0.5 + 1e-4 * rnorm(n_walks)        
        }
        for (i in 1:3)
        {
            # Loop over scale parameters of bias model for each output
            post_samples[(num_model_params+i),,1] = 
                init_sd_bias[[trc_outputs[i]]] + 
                sd_multiplier[[trc_outputs[i]]] * rnorm(n_walks)
        }        
    }

    return(post_samples)
}
