#
# title     : RunTrcMCMC.R
# purpose   : R script with run MCMC AIES sampler on the FEBA model
# author    : WD41, LRS/EPFL/PSI
# date      : July 2017
#
# Load Library ----------------------------------------------------------------
library(DiceKriging)
library("optparse")
library("rgw")

# Source Auxiliary Functions
source("./r-scripts/sourceAuxFunctions.R")
# Source Helper Functions
source("./r-scripts/sourceHelpers.R")

# Make and parse command line arguments ---------------------------------------
option_list <- list(
    make_option(c("--trc_data_bias"), type = "character", default = NULL,
                help = "Bias simulation campaign RDS",
                metavar = "character"),
    make_option(c("--pca"), type = "character", default = NULL,
                help = "PCA results RDS",
                metavar = "character"),
    make_option(c("--pc_gp"), type = "character", default = NULL,
                help = "GP metamodels for PC Scores RDS",
                metavar = "character"),
    make_option(c("--bias_gp"), type = "character", default = NULL,
                help = "GP model for bias RDS",
                metavar = "character"),
    make_option(c("--output"), type = "character", default = NULL,
                help = "The output (tc | dp | co | all)",
                metavar = "character"),
    make_option(c("-n", "--numprocs"), type = "integer", default = 1,
                help = "The number of processors",
                metavar = "integer"),
    make_option(c("--n_walks"), type = "integer", default = 20,
                help = "The number of walkers",
                metavar = "integer"),
    make_option(c("--n_iters"), type = "integer", default = 10,
                help = "The number of iterations",
                metavar = "integer")
)

opt_parser <- OptionParser(option_list = option_list) # Create opt. parser obj.
opt <- parse_args(opt_parser)                         # Parse the object

pca_files <- strsplit(opt$pca, ",")[[1]]
gp_bias_files <- strsplit(opt$bias_gp, ",")[[1]]
gps_pcs_files <- strsplit(opt$pc_gp, ",")[[1]]

# Read Data -------------------------------------------------------------------
if (opt$output == "all")
{
    trc_outputs <- c("tc", "dp", "co")
} else
{
    trc_outputs <- opt$output
}

trc_data_bias <- readRDS(opt$trc_data_bias) # Compiled TRACE runs

trc_pca <- list()
trc_gp_bias <- list()
trc_gps_pcs <- list()
for (i in 1:length(trc_outputs))
{
    trc_pca[[trc_outputs[i]]] <- readRDS(trimws(pca_files[i]))
    trc_gp_bias[[trc_outputs[i]]] <- readRDS(trimws(gp_bias_files[i]))
    trc_gps_pcs[[trc_outputs[i]]] <- readRDS(trimws(gps_pcs_files[i]))
}

# Global variables ------------------------------------------------------------
num_pc <- list("tc" = 7, "dp" = 10, "co" = 5)

# Pre-processing --------------------------------------------------------------
init_sd_bias <- list()
time_idx <- list()
exp_vec <- list()
trc_pca_ave <- list()
trc_pca_lds <- list()
for (trc_output in trc_outputs)
{
    # Get the initial estimate of scale parameter of the bias model
    init_sd_bias[[trc_output]] <- sqrt(
        trc_gp_bias[[trc_output]]@covariance@sd2)
    # Mean output and PC loadings
    trc_pca_ave[[trc_output]] <- trc_pca[[trc_output]]$yy_train_ave
    trc_pca_lds[[trc_output]] <- trc_pca[[trc_output]]$yy_train_loadings  
    
    if (trc_output == "tc")
    {
        # Clad temperature output, restricted time
        time_idx[["tc"]] <- GetTimeIdxTC(trc_data_bias)
        # Clad temperature, experimental data
        exp_vec[["tc"]] <- GetExpDataTC(trc_data_bias)
    } else if (trc_output == "dp")
    {
        # Pressure Drop output, restricted time
        time_idx[["dp"]] <- GetTimeIdxDP(trc_data_bias)
        # Pressure drop, experimental data
        exp_vec[["dp"]] <- GetExpDataDP(trc_data_bias)
    } else if (trc_output == "co")
    {
        # Liquid carryover output, restricted time
        time_idx[["co"]] <- GetTimeIdxCO(trc_data_bias)
        # Liquid carryover, experimental data
        exp_vec[["co"]] <- GetExpDataCO(trc_data_bias)
    }    
}

# Setup the probabilistic model -----------------------------------------------
# Log Likelihood
if (opt$output == "all")
{
    exp_idx_max_tc <- GetExpTimeQuenchIdx(trc_data_bias)
    exp_idx_max_dp <- rep(length(trc_data_bias$exp_data[[2]][,1]), 4)
    
    # Now the input parameters becomes 15 (12 + 3 bias scale parameters)
    llk <- function(x) {
        GetLogLikelihood(c(x[1:12], x[13]), 
                      exp_vec[["tc"]], num_pc[["tc"]], time_idx[["tc"]], 
                      trc_gps_pcs[["tc"]], trc_pca_ave[["tc"]], 
                      trc_pca_lds[["tc"]], trc_gp_bias[["tc"]], 
                      ax_locs = unique(trc_gp_bias[["tc"]]@X[,5]),
                      time_pts = unique(trc_gp_bias[["tc"]]@X[,6])[-1],
                      exp_idx_max_tc - 1) +
            GetLogLikelihood(c(x[1:12], x[14]), 
                          exp_vec[["dp"]], num_pc[["dp"]], time_idx[["dp"]], 
                          trc_gps_pcs[["dp"]], trc_pca_ave[["dp"]], 
                          trc_pca_lds[["dp"]], trc_gp_bias[["dp"]], 
                          ax_locs = unique(trc_gp_bias[["dp"]]@X[,5]), 
                          time_pts = unique(trc_gp_bias[["dp"]]@X[,6]), 
                          exp_idx_max_dp) + 
            GetLogLikelihood(c(x[1:12], x[15]), 
                          exp_vec[["co"]], num_pc[["co"]], time_idx[["co"]], 
                          trc_gps_pcs[["co"]], trc_pca_ave[["co"]], 
                          trc_pca_lds[["co"]], trc_gp_bias[["co"]],
                          ax_locs = NULL, 
                          time_pts = unique(trc_gp_bias[["co"]]@X[,5]),
                          exp_idx_max = NULL)
    }
} else if (trc_outputs[1] == "tc")
{
    ax_locs <- unique(trc_gp_bias[["tc"]]@X[,5])
    time_pts <- unique(trc_gp_bias[["tc"]]@X[,6])
    exp_idx_max <- GetExpTimeQuenchIdx(trc_data_bias)
    llk <- function(x) {
        GetLogLikelihood(x, 
                      exp_vec[["tc"]], num_pc[["tc"]], time_idx[["tc"]], 
                      trc_gps_pcs[["tc"]], trc_pca_ave[["tc"]], 
                      trc_pca_lds[["tc"]], trc_gp_bias[["tc"]], 
                      ax_locs, time_pts[-1], exp_idx_max - 1)
    }
} else if (trc_outputs[1] == "dp")
{
    ax_locs <- unique(trc_gp_bias[["dp"]]@X[,5])
    time_pts <- unique(trc_gp_bias[["dp"]]@X[,6])
    exp_idx_max <- rep(length(trc_data_bias$exp_data[[2]][,1]), 4)
    
    llk <- function(x) {
        GetLogLikelihood(x, 
                      exp_vec[["dp"]], num_pc[["dp"]], time_idx[["dp"]], 
                      trc_gps_pcs[["dp"]], trc_pca_ave[["dp"]], 
                      trc_pca_lds[["dp"]], trc_gp_bias[["dp"]], 
                      ax_locs, time_pts, exp_idx_max)
    }
} else if (trc_outputs[1] == "co")
{
    ax_locs <- NULL
    time_pts <- unique(trc_gp_bias[["co"]]@X[,5])

    llk <- function(x) {
        GetLogLikelihood(x, 
                      exp_vec[["co"]], num_pc[["co"]], time_idx[["co"]], 
                      trc_gps_pcs[["co"]], trc_pca_ave[["co"]], 
                      trc_pca_lds[["co"]], trc_gp_bias[["co"]],
                      ax_locs, time_pts)
    }
}

# Log Prior
# First 12, Model parameters -> uniform [0, 1]
# Bias scale parameters (standard deviation) -> estimate from the GP model
if (opt$output != "all")
{
    lprior <- function(x)
    {
        return(sum(dunif(x[1:12], 0, 1, log = T)) + 
                   dhalfcauchy(x[13], init_sd_bias[[1]], log = T))
    }
} else
{
    lprior <- function(x)
    {
        return(sum(dunif(x[1:12], 0, 1, log = T)) + 
                   dhalfcauchy(x[13], init_sd_bias[[1]], log = T) +
                   dhalfcauchy(x[14], init_sd_bias[[2]], log = T) +
                   dhalfcauchy(x[15], init_sd_bias[[3]], log = T))
    }
}

# Log Posterior
lpost <- function(x)
{
    lp <- lprior(x)
    if (is.infinite(lp))
    {
        return(-Inf)
    } else
    {
        return(lp + llk(x))
    }
}

# MCMC Sampler ----------------------------------------------------------------
# Set up AIES Walkers
n_walks <- opt$n_walks
n_iters <- opt$n_iters
sd_multiplier <- list("tc" = 10., "dp" = 100., "co" = 0.1)

# Set up initial values
if (opt$output != "all")
{
    post_samples <- array(NA, dim = c(13, n_walks, n_iters))
    # Set up initial values, 13 parameters, single output
    for (i in 1:12)
    {
        post_samples[i,,1] = 0.5 + 1e-4 * rnorm(n_walks)
        
    }
    post_samples[13,,1] = init_sd_bias[[1]] + 
        sd_multiplier[[trc_outputs[1]]] * rnorm(n_walks)
} else
{
    post_samples <- array(NA, dim = c(15, n_walks, n_iters))
    # Set up initial values, 15 parameters, multiple outputs
    for (i in 1:12)
    {
        post_samples[i,,1] = 0.5 + 1e-4 * rnorm(n_walks)
        
    }
    post_samples[13,,1] = init_sd_bias[[1]] + 
        sd_multiplier[[trc_outputs[1]]] * rnorm(n_walks)
    post_samples[14,,1] = init_sd_bias[[2]] + 
        sd_multiplier[[trc_outputs[2]]] * rnorm(n_walks)
    post_samples[15,,1] = init_sd_bias[[3]] + 
        sd_multiplier[[trc_outputs[3]]] * rnorm(n_walks)
}

# Run Sampler
ens <- GoodmanWeare.rem(post_samples, lpost, 
    mc.cores = opt$numprocs, mention.every = 10)

# Save the samples
saveRDS(ens, 
  paste("ens-", opt$output, "-", n_walks, "-", n_iters, ".Rds", sep = ""))
