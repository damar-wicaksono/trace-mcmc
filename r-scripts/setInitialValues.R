#' Set initial values of the chain in the walker
#'
#' @param num_model_params Numeric. The number of model parameters
#' @param fix_scale Logical. Flag whether to fix scale
#' @param init_sd_bias List. The iniitial bias model scale parameter(s)
#' @param n_walks Numeric. Number of walkers
#' @param n_iters Numeric. Number of iterations
#' @param trc_outputs Character vector. The type of outputs to be considered
#' @return 3D array, 1st dimension is the container for each parameter, 
#'  2nd dimension is the container for each walker,
#'  3rd dimension is the sample per each iteration
setInitialValues <- function(num_model_params, fix_scale,
                             init_sd_bias, n_walks, n_iters)
{
    sd_multiplier <- list("tc" = 10., "dp" = 100., "co" = 0.1)

    if (fix_scale)
    {
        num_params <- num_model_params
    } else
    {
        num_params <- num_model_params + length(init_sd_bias)
    }

    post_samples <- array(NA, dim = c(num_params, n_walks, n_iters))
    
    # Set up initial values for model parameters
    for (i in 1:num_model_params)
    {
        # Initial values at nominal then perturbed
        post_samples[i,,1] = 0.5 + 1e-3 * rnorm(n_walks)
    }

    # Set up initial values for scale parameters
    if (!fix_scale)
    {   
        for (i in 1:length(init_sd_bias))
        {
            # Loop over scale parameters of bias model for each output
            post_samples[(num_model_params+i),,1] = 
                init_sd_bias[[i]] + 
                sd_multiplier[[names(init_sd_bias[i])]] * rnorm(n_walks)
        }
    }

    return(post_samples)
}
