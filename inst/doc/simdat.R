## ---- echo=TRUE, message=FALSE, warning=FALSE, cache=FALSE---------------
library(devtools)
install_github("TobiasRoth/detecfilter")
library(detecfilter)

## ---- cache=TRUE---------------------------------------------------------
# Calculate observed community matrix
z_obs <- apply(plantsBDM$y, c(1,2), max)

# Calculate observed occupancies for all species
psi <- apply(z_obs, 2, mean)

# Mean and sd at logit-scale
(mean.lpsi <- mean(qlogis(psi)))
(sig.lpsi <- sd(qlogis(psi)))

## ---- fig.height=4, cache=TRUE-------------------------------------------
par(mfrow = c(1,2))
hist(psi, breaks = 100, ylim = c(0,500), main = "", 
     ylab = "Number of species", xlab = "Occupancy (observed)")
hist(plogis(rnorm(1733, mean.lpsi, sig.lpsi)), 
     breaks = 100, ylim = c(0,500), main = "",
     ylab = "Number of species", xlab = "Occupancy (simulated)")

## ---- fig.height=4, cache=TRUE-------------------------------------------
# Choose how many species should be simulated
nspec <- 50  

# Choose the values for the four parameters to describe species' occupancies along gradient
mean.lpsi <- 1.0        # Mean intercept of species' occupancies
sigma.lpsi <- 1.0       # SD of species differences in intercept
mu.beta.lpsi <- -0.3    # Mean slope of species' occupancies along gradient
sig.beta.lpsi  <- 0.5   # SD of species differences in slopes

# Simulate intercept (b0) and slope (b1) for each species
b0 <- rnorm(nspec, mean.lpsi, sigma.lpsi) 
b1 <- rnorm(nspec, mu.beta.lpsi, sig.beta.lpsi)

# Plot species occupancies along gradient
plot(NA, xlim = c(-3, 3), ylim = c(0, 1), 
     xlab = "gradient (e.g. elevation)", ylab = "Occupancy")
gradient <- seq(-3, 3, 0.1)
for(k in 1:nspec) {
  lpsi <- b0[k] + b1[k] * gradient 
  points(gradient, plogis(lpsi), lwd = 0.5, ty = "l", col = "grey")
}

## ---- cache=TRUE---------------------------------------------------------
# Number of sites in the meta-community 
nsite <- 100

# Distribution of sites along the gradient
# Note that the choice of the distribution is not of particular importance, and 
# we may also use a uniform distiruntion
gradient <- rnorm(nsite, 0, 1.5) 

# Simulation of the realized meta-community
z <- matrix(NA, nsite, nspec)
for(k in 1:nspec) {
  lpsi <- b0[k] + b1[k] * gradient 
  z[,k] <- rbinom(nspec, 1, plogis(lpsi))
}

## ---- cache=TRUE---------------------------------------------------------
# Distribution of functional trait expresssion across species
# (we assume that we standardized the trait values and the range 
# of trait values covers about 6 standard deviations)
trait <- rnorm(nspec, 0, 1.5)

# Calculate community mean
CM <- apply(z, 1, function(x) mean(trait[x==1]))

# Plot CMs along gradient
plot(gradient, CM, pch = 16, cex = 0.7, ylim = c(-1.5, 1.5))
abline(lm(CM ~ gradient))
abline(h = mean(trait), lty = 2)

## ------------------------------------------------------------------------
# Effect of functional trait on species response to gradient
mu.FTfilter.lpsi <- -1

# Species response (slope) also depends on functional trait
b1 <- rnorm(nspec, mu.beta.lpsi + mu.FTfilter.lpsi * trait, sig.beta.lpsi) 

# Simulation of the realized meta-community (same code as above)
z <- matrix(NA, nsite, nspec)
for(k in 1:nspec) {
  lpsi <- b0[k] + b1[k] * gradient 
  z[,k] <- rbinom(nspec, 1, plogis(lpsi))
}

# Calculate community mean
CM <- apply(z, 1, function(x) mean(trait[x==1]))

# Plot CMs along gradient
plot(gradient, CM, pch = 16, cex = 0.7, ylim = c(-1.5, 1.5),
     xlab = "Gradient (e.g. elevation)", ylab = "CM (e.g. average plant height)")
abline(lm(CM~gradient))
abline(h = mean(trait), lty = 2)

## ---- fig.height=4, cache=TRUE-------------------------------------------
# Choose the values for the four parameters to describe species' detection probability
mean.lp <- 0.8        # Mean intercept of species' detection probability
sigma.lp <- 1.0       # SD of species differences in intercept
mu.beta.lp <- 1       # Mean species' phenological effect on detection probability
sig.beta.lp  <- 1     # SD of species differences in phenological effect

# Simulate intercept (b0) and slope (b1) for each species
a0 <- rnorm(nspec, mean.lp, sigma.lp) 
a1 <- rnorm(nspec, mu.beta.lp, sig.beta.lp)

# Distribution of data of surveys
date <- round(rnorm(nspec, 0, 10), 0)

# Simulation of observed community
y <- matrix(NA, nsite, nspec)
for(k in 1:nspec) {
  lp <- a0[k] + a1[k] * date 
  y[,k] <- rbinom(nspec, 1, z[,k]*plogis(lpsi))
}

## ---- fig.height=4-------------------------------------------------------
# Calculate species richness
SR_true <- apply(z, 1, sum)
SR_obs <- apply(y, 1, sum)

# Calculate community mean of trait expression
CM_true <- apply(z, 1, function(x) mean(trait[x==1]))
CM_obs <- apply(y, 1, function(x) mean(trait[x==1]))

# Make plots
par(mfrow = c(1,2))
plot(gradient, SR_obs, pch = 16, cex = 0.7, ylab = "Species richness", ylim = c(0,50))
abline(lm(SR_obs ~ gradient))
points(gradient, SR_true, cex = 0.7, col = "orange")
abline(lm(SR_true ~ gradient), col = "orange")
plot(gradient, CM_obs, pch = 16, cex = 0.7, ylab = "CM", ylim = c(-2,2))
abline(lm(CM_obs ~ gradient))
points(gradient, CM_true, cex = 0.7, col = "orange")
abline(lm(CM_true ~ gradient), col = "orange")


## ---- fig.height=4, cache=TRUE-------------------------------------------
# Effect of functional trait on the average detection probability of a species
mu.FTfilter.lp <- 2

# Average detection probability (intercept) also depends on functional trait 
a0 <- rnorm(nspec, mean.lp + mu.FTfilter.lp * trait,  sigma.lp) 

# Simulation of observed community (same code as above)
y <- matrix(NA, nsite, nspec)
for(k in 1:nspec) {
  lp <- a0[k] + a1[k] * date 
  y[,k] <- rbinom(nspec, 1, z[,k]*plogis(lp))
}

# Calculate community mean of trait expression
CM_obs <- apply(y, 1, function(x) mean(trait[x==1]))

# Make plots
par(mfrow = c(1,2))
plot(trait, plogis(a0), pch = 16, cex = 0.7, ylim = c(0,1),
     ylab = "Avg.  detection probability", xlab = "Trait (e.g. plant height)")
plot(gradient, CM_obs, pch = 16, cex = 0.7, ylab = "CM", ylim = c(-2,2))
abline(lm(CM_obs ~ gradient))
points(gradient, CM_true, cex = 0.7, col = "orange")
abline(lm(CM_true ~ gradient), col = "orange")


