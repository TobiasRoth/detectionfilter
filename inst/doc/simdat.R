## ---- echo=TRUE, message=FALSE, warning=FALSE, cache=FALSE---------------
library(devtools)
install_github("TobiasRoth/detectionfilter")
library(detectionfilter)

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


