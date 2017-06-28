## ---- echo=TRUE, message=FALSE, warning=FALSE, cache=TRUE----------------
library(devtools)
install_github("TobiasRoth/detectionfilter")
library(detectionfilter)

## ---- cache=TRUE---------------------------------------------------------
set.seed(1234)
dat <- simcom(mu.FTfilter.lpsi = -0.5, mean.psi = 0.8, 
              mu.FTfilter.lp = 1, mean.p = 0.8,
              nsite = 200, nspec = 100, nrep = 2)

## ---- cache=TRUE---------------------------------------------------------
# Calculate observed meta-community (i.e. merge the observations from the two visits)
commat_obs <- apply(dat$y, c(1,2), max)

# Show species with few observations
sort(apply(commat_obs, 2, sum))[1:10]

## ---- message=FALSE, warning=FALSE, cache=TRUE---------------------------
# Package to fit hierarchical models of occurrence and abundance data
library(unmarked)

# Prepare detection-corrected meta-community matrix
z <- array(NA, dim = c(dim(dat$y)[1], dim(dat$y)[2]))

# Apply hierarchical model to all species seperate
for(k in 1:dim(dat$y)[2]) {
  d <- unmarked::unmarkedFrameOccu(y = dat$y[,k,], obsCovs = list(date = dat$date),
                         siteCovs = data.frame(gradient = dat$gradient))
  res <- unmarked::occu(~ date ~ gradient, data = d, se = TRUE)
  z[,k] <- unmarked::bup(unmarked::ranef(res), stat = "mode")
}

## ---- cache=TRUE, fig.height=5-------------------------------------------
# Function to calculate mean trait expression of species in a community (CM)
CM <- function(x) mean(dat$traitmat[x])

# Calculate CMs from observed meta-community matrix, from detection-corrected 
# meta-community matrix and from true mata-community matrix
CM_obs <- apply(commat_obs==1, 1, CM)
CM_cor <- apply(z>0.5, 1, CM)
CM_true <- apply(dat$z_true==1, 1, CM)

# Plot CM along gradient
plot(dat$gradient, CM_obs, pch = 16, cex = 0.7, ylim = c(-1, 1), 
     xlab = "Gradient", ylab = "Community mean")
points(dat$gradient, CM_true, cex = 0.7, col = "orange")
points(dat$gradient, CM_cor, pch = "+", cex = 0.7, col = "red")
abline(lm(CM_obs ~ dat$gradient))
abline(lm(CM_true ~ dat$gradient), col = "orange")
abline(lm(CM_cor ~ dat$gradient), col = "red")

## ---- cache=TRUE---------------------------------------------------------
set.seed(1234)
dat <- simcom(mu.FTfilter.lpsi = -0.5, mean.psi = 0.5, 
              mu.FTfilter.lp = 2, mean.p = 0.5,
              nsite = 200, nspec = 100, nrep = 2)

# Observed community matrix
commat_obs <- apply(dat$y, c(1,2), max)

# Number of observed species
ncol(commat_obs)

## ---- cache=TRUE---------------------------------------------------------
sum(apply(dat$z_true, 2, sum)>0)

## ------------------------------------------------------------------------
# Proportion of  species observed at least once but in less than 10% of sites
mean(apply(commat_obs, 2, sum) < 20)

## ---- message=FALSE, warning=FALSE, cache=TRUE, fig.height=5-------------
# Estimate detection-corrected meta-community matrix
z <- array(NA, dim = c(dim(dat$y)[1], dim(dat$y)[2]))
for(k in 1:dim(dat$y)[2]) {
  d <- unmarked::unmarkedFrameOccu(y = dat$y[,k,], obsCovs = list(date = dat$date),
                         siteCovs = data.frame(gradient = dat$gradient))
  res <- unmarked::occu(~ date ~ gradient, data = d, se = TRUE)
  z[,k] <- unmarked::bup(unmarked::ranef(res), stat = "mode")
}

# Calculate CMs
CM_obs <- apply(commat_obs==1, 1, CM)
CM_cor <- apply(z>0.5, 1, CM)
CM_true <- apply(dat$z_true==1, 1, CM)

# Plot CMs along gradient
plot(dat$gradient, CM_obs, pch = 16, cex = 0.7, ylim = c(-1, 1), 
     xlab = "Gradient", ylab = "Community mean")
points(dat$gradient, CM_true, cex = 0.7, col = "orange")
points(dat$gradient, CM_cor, pch = "+", cex = 0.7, col = "red")
abline(lm(CM_obs ~ dat$gradient))
abline(lm(CM_true ~ dat$gradient), col = "orange")
abline(lm(CM_cor ~ dat$gradient), col = "red")

