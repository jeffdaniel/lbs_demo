# ------------------------------------------------------------------------------
#             Efficient modelling of presence-only species data via
#                           local background sampling
# ------------------------------------------------------------------------------
#
# This script demonstrates local background sampling with simulated data. 
# The procedure can be easily implemented in R using the 'spatstat' package,
# with an important caveat explained below.
#
# ------------------------------------------------------------------------------

# Load 'spatstat' package.
library(spatstat)

# Load simulated data, consisting of two objects:
#    - X, a spatial point pattern,
#    - data, a list of five pixel images serving as spatial covariates.
load("data.Rda")

# Create a binary image mask covering the spatial domain.
mask <- as.mask(data[[1]])

# Set total number of background points.
m <- 5000

# Fit a pilot model via conditional logistic regression.
S_pilot   <- runifpoint(m, mask)
Q_pilot   <- quadscheme.logi(X, S_pilot)
fit_pilot <- ppm(Q_pilot ~ ., data = data, method = "logi")

# Obtain pilot intensity estimate.
lambda_pilot <- predict(fit_pilot, locations = mask)

# Obtain offset term.
offset_pilot <- m * lambda_pilot / integral.im(lambda_pilot)
      
# Simulate local background point pattern.
S_local <- rpoint(m, lambda_pilot)
Q_local <- quadscheme.logi(X, S_local)

# NOTE: When used to fit models via conditional logistic regression, the 'ppm'
# function assumes that the background point process has homogeneous intensity
# and includes the offset -log(rho), where rho is the total number of background
# points divided by the area of the spatial domain. To correctly implement the
# local background sampling procedure, we must remove this offset term. Here we
# do so by setting the value of rho to 1, yielding an offset of 0.
Q_local$param$rho <- 1

# Fit final model via conditional logistic regression with offset.
fit_local <- ppm(Q_local ~ . + offset(-log(offset_pilot)), data = data, 
                 method = "logi")

# NOTE: To obtain correct intensity estimates, the result computed by
# 'predict.ppm' must be multiplied by the offset term.
lambda_final <- predict(fit_local, locations = mask) * offset_pilot
plot(log(lambda_final))
