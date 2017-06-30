# RGW Package
install.packages("rgw")
library(rgw)

# In this example, we'll sample from a simple 2D Gaussian.

# Define the log-posterior function
lnP = function(x) sum( dnorm(x, c(0,1), c(pi, exp(0.5)), log=TRUE) )

# Initialize an ensemble of 100 walkers. We'll take 100 steps, saving the ensemble after each.
nwalk = 100
post = array(NA, dim=c(2, nwalk, 1001))
post[1,,1] = rnorm(nwalk, 0, 0.1)
post[2,,1] = rnorm(nwalk, 1, 0.1)

# Run
post = GoodmanWeare.rem(post, lnP)

# Plot the final ensemble
plot(post[1,,101], post[2,,101])
# Look at the trace of each parameter for one of the walkers.
plot(post[1,1,], type = "l", col = rgb(0, 0, 0, 0.25))
for (i in 2:nwalk) lines(post[1,i,], type = "l", col = rgb(0, 0, 0, 0.25))
abline(h = 0, col = "blue", lwd = 2)
plot(post[2,1,], type = "l", col = rgb(0, 0, 0, 0.25))
for (i in 2:nwalk) lines(post[2,i,], type = "l", col = rgb(0, 0, 0, 0.25))
abline(h = 1, col = "blue", lwd = 2)

# Go on to get confidence intervals, make niftier plots, etc.
dim(post)
plot(post[1,,101], type = "l")
hist(post[1,,101])
summary(post[1,,101])

abline(v = 0)
hist(post[2,,101])
abline(v = 1)
plot(post[2,1,], type = "l")

# ---
# In this example, we'll sample from a simple 2D Gaussian
# Define the log-posterior function
lnP = function(x) sum( dnorm(x, c(0,1), c(pi, exp(0.5)), log=TRUE) )
# Initialize an ensemble of 100 walkers
nwalk = 10000
ensemble = array(dim=c(2, nwalk))
ensemble[1,] = rnorm(nwalk, 0, 0.1)
ensemble[2,] = rnorm(nwalk, 1, 0.1)
# Run for a bit
ens2 = GoodmanWeare(ensemble, lnP, 100, mc.cores=1)
# Plot the resulting ensemble
plot(t(ens2$ensemble))
# Compare to a direct draw from the posterior distribution
points(rnorm(nwalk, 0, pi), rnorm(nwalk, 1, exp(0.5)), col=2, pch=3)

# Fitting Data Example --------------------------------------------------------
# Originally appeared as Python example in the emcee main website
set.sed(666)

# Setup the true parameter
m_true <- -0.9594
b_true <- 4.294
f_true <- 0.534

# Generate synthetic data
N <- 50                                 # The total number of observation
x <- sort(10 * runif(N))                # The observation point
y_err <- 0.1 + 0.5 * runif(N)           # The true observation error
y <- m_true * x + b_true                # The true model
y <- y + abs(f_true * y) * rnorm(N)    # Introduce a bias to the model
y <- y + y_err * rnorm(N)               # Further distort the true error

plot(x, y, ylim = c(-6, 6))
lines(x, m_true * x + b_true)
for (i in 1:length(y_err))
{
    points(x[i], y[i] + y_err[i], pch = "-")
    points(x[i], y[i] - y_err[i], pch = "-")
}

# Least Square solution
A <- cbind(rep(1, length(x)), x)
C <- diag(y_err * y_err)
covar <- solve(t(A) %*% solve(C) %*% A)
bm <- covar %*% t(A) %*% solve(C) %*% y
b_ls <- bm[1]
m_ls <- bm[2]
abline(a = b_ls, b = m_ls, lwd = 2, lty = 2)

# Setup Log-Prior
lnprior <- function(theta)
{
    m <- theta[1]
    b <- theta[2]
    #lnf <- theta[3]
    
    return(dunif(m, -5, 0.5, log = T) + dunif(b, 0, 10, log = T)) 
#           dunif(lnf, -10, 1.0, log = T))
    #if ((m > -5 && m < 0.5) && (b > 0.0 && b < 10.0) && (lnf > -10.0 && lnf < 1.0))
    #{
    #    return(0.0)
    #} else
    #{
    #    return(-Inf)
    #}
}

# Setup Log-likelihood
, x = x, y = y, yerr = y_err
lnlike <- function(theta)
{
    m <- theta[1]
    b <- theta[2]
    lnf <- theta[3]
    
    model <- m * x + b
    inv_sigma2 = 1.0 / (y_err^2 + model^2 * exp(2 * lnf))
    
    return(-1 * -0.5 * sum((y - model)^2 * inv_sigma2 - log(inv_sigma2)))
}

# Setup Log-posterior
lnpost <- function(theta)
{
    lp <- lnprior(theta)
    if (is.infinite(lp))
    {
        return(-Inf)
    } else
    {
        return(lp + lnlike(theta, x, y, y_err))
    }
}

optim(par = c(-0.9, 4, -1), lnlike, method = "BFGS")

x <- c(0.17006809,  0.31992597,  0.79839803,  0.83969083,  0.85521133,
0.86430831,  0.93204473,  1.15041274,  1.30420486,  1.34773303,
1.36125441,  1.45923308,  1.56367109,  2.18749753,  2.42498413,
2.60606778,  3.19119707,  3.63817337,  3.75035431,  4.2381022 ,
4.29814407,  4.56346462,  4.59200088,  4.78770972,  4.93046958,
5.12921474,  5.35396938,  5.56414315,  5.73025443,  6.40152157,
6.57345799,  6.63758379,  6.65233434,  6.73960383,  7.14855074,
7.73546755,  7.79315669,  7.8683726 ,  8.0027163 ,  8.00528833,
8.0249634 ,  8.06892082,  8.20307725,  8.84191015,  8.96601818,
9.32032505,  9.35686615,  9.52783263,  9.73922635,  9.94463153)

y <- c(7.73772992,  5.49693267,  4.39095111,  3.32618534,  1.54645103,
       1.91597314,  2.10663586,  3.61641386,  3.791394  ,  1.86140523,
       2.8794792 ,  4.38553385, -0.83969421,  3.5049713 ,  2.73121819,
       2.8814837 ,  0.89496683,  0.78846637, -0.08938221,  0.22745595,
       0.28445952,  0.03050963, -0.1256389 , -0.28884368, -0.32011869,
       0.08340355, -2.10550927, -1.18363358, -2.28727491, -3.23159607,
       -2.46765836, -2.61414535, -0.17753437, -2.72467123, -0.31249129,
       -4.03804459, -1.37334798, -3.70076924, -3.33651723, -4.432594  ,
       -4.43680834, -5.03852925, -2.03530896, -1.96572261, -4.90282249,
       -7.19910215, -6.17213199, -3.76840964, -5.28674698, -4.72181485)
y_err <- c(
    0.58633756,  0.45824509,  0.54513126,  0.23920054,  0.52989702,
    0.14921012,  0.58738922,  0.15791897,  0.27621634,  0.48926804,
    0.29960707,  0.1549446 ,  0.30406605,  0.38292939,  0.29643064,
    0.54686075,  0.27946605,  0.11215477,  0.24107133,  0.42463668,
    0.21312195,  0.23090935,  0.23728126,  0.27946674,  0.13784194,
    0.47973564,  0.50287616,  0.19968825,  0.44093171,  0.50622733,
    0.5871678 ,  0.59143584,  0.47833901,  0.20457842,  0.13633574,
    0.1419815 ,  0.10706326,  0.58360097,  0.40583197,  0.20476346,
    0.47523769,  0.50110128,  0.44290436,  0.57150628,  0.4289756 ,
    0.21551228,  0.11486808,  0.49509282,  0.52559147,  0.10527325
)
# Setup AIES Walker
nwalk = 12
post_samples = array(NA, dim=c(2, nwalk, 5000))
post_samples[1,,1] = rep(0.1, nwalk) #runif(nwalk, -5, 0.5)
post_samples[2,,1] = rep(0.1, nwalk) #runif(nwalk, 0, 10.0)
post_samples[3,,1] = rep(, nwalk) #runif(nwalk, -10., 1)

# Run
ens3 = GoodmanWeare.rem(post_samples, lnpost, mc.cores = 1, mention.every = 1000)

res <- matrix(ens3[,,1000:5000], nrow = 2)

# parameter m, trace plot
plot(ens3[1,1,], type = "l", ylim = c(-5, 0.5))
for (i in 2:nwalk) lines(ens3[1,i,], type = "l")
abline(h = m_true, col = "red")
# parameter b, trace plot
plot(ens3[2,1,], type = "l", ylim = c(0, 10))
for (i in 2:nwalk) lines(ens3[2,i,], type = "l")
abline(h = b_true, col = "red")
# parameter log(f), trace plot
plot(ens3[3,1,], type = "l", ylim = c(-10, 1))
for (i in 2:nwalk) lines(ens3[3,i,], type = "l")
abline(h = log(f_true), col = "red")

plot(as.vector(ens3[1,,50:101]), type = "l")
plot(res[2,2], type = "l")
abline(h = b_true, col = "red")
plot(as.vector(ens3[3,,50:101]), type = "l")
abline(h = log(f_true), col = "red")

hist(res[1,], breaks = seq(-1,-0.5,length.out = 100))
abline(v = m_true, col = "red")
hist(res[2,], breaks = seq(0,10,length.out = 100))
abline(v = b_true, col = "red")
hist(res[3,], breaks = seq(-10,1,length.out = 100))
abline(v = log(f_true), col = "red")

plot(res[1,], res[2,], col = rgb(0,0,0,0.01))
z <- kde2d(res[1,], res[2,], n = 50)
contour(z, drawlabels = F, nlevels = 11, add = T, col = "red")
plot(res[1,], res[3,], col = rgb(0,0,0,0.01))
plot(res[2,], res[3,], col = rgb(0,0,0,0.01), ylim = c(-10, 1), xlim = c(0, 10))
z <- kde2d(res[2,], res[3,], n = 50)
contour(z, drawlabels = F, nlevels = 11, add = T, col = "red")

contour(z, drawlabels = F, nlevels = 11, add = T, col = "red")
smoothScatter(res[1:2,])
plot(1:60, ens3$ensemble[1,], type = "l")
z <- kde2d(res[1,], res[2,], n = 50)

, x, y, yerr