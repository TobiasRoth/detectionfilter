## ---- echo=TRUE, message=FALSE, warning=FALSE, cache=FALSE---------------
library(devtools)
install_github("TobiasRoth/detecfilter")
library(detecfilter)

## ---- message=FALSE, warning=FALSE, cache=TRUE---------------------------
library(unmarked)
z <- array(NA, dim = c(dim(dat$y)[1], dim(dat$y)[2]))
a0 <- numeric(dim(dat$y)[2])
a0[!is.na(a0)] <- NA
b1 <- numeric(dim(dat$y)[2])
b1[!is.na(b1)] <- NA
for(k in 1:dim(dat$y)[2]) {
  d <- unmarkedFrameOccu(y = dat$y[,k,], siteCovs = data.frame(gradient = dat$gradient))
  try({
    res <- occu(~ 1 ~ gradient, data = d, se = TRUE)
    a0[k] <- plogis(coef(res)["p(Int)"])
    b1[k] <- coef(res)["psi(gradient)"]
    z[,k] <- bup(ranef(res), stat = "mode")
 })
}

## ------------------------------------------------------------------------
# Calculate community mean value of trait expression (CM)
commat_obs <- apply(dat$y, c(1,2), max)
CM_obs <- apply(commat_obs, 1, function(x) mean(dat$traitmat[dat$detected.at.all][x==1]))
CM_cor <- apply(z>0.5, 1, function(x) mean(dat$traitmat[dat$detected.at.all][x], na.rm = TRUE))
CM_true <- apply(dat$z_true, 1, function(x) mean(dat$traitmat[x==1]))

# Plot CM along gradient
plot(dat$gradient, CM_obs, pch = 16, cex = 0.7, ylim = c(-1.5, 1.5))
points(dat$gradient, CM_true, cex = 0.7, col = "orange")
points(dat$gradient, CM_cor, pch = "+", cex = 0.7, col = "red")
abline(lm(CM_obs ~ dat$gradient))
abline(lm(CM_true ~ dat$gradient), col = "orange")
abline(lm(CM_cor ~ dat$gradient), col = "red")
abline(h=mean(dat$traitmat), lty = 2)
abline(h=0)

## ---- eval=FALSE, include=FALSE------------------------------------------
#  plot(dat$P[dat$detected.at.all], a0)
#  abline(0,1)
#  plot(dat$traitmat[dat$detected.at.all], a0)
#  
#  plot(dat$traitmat[dat$detected.at.all], b1, ylim = c(-2,2))
#  abline(h=0)
#  
#  plot(dat$gradient, apply(commat_obs, 1, sum), pch = 16, cex = 0.7, ylim = c(0, 40))
#  points(dat$gradient, apply(dat$z_true, 1, sum), cex = 0.7, col = "orange")
#  points(dat$gradient, apply(z>0.5, 1, sum), pch = "+", cex = 0.7, col = "red")
#  abline(h=0)
#  abline(h=mean(dat$traitmat), lty = 2)
#  
#  # Calculate community mean value of trait expression (CM)
#  commat_obs <- apply(dat$y, c(1,2), max)
#  ausw <- b1 < 10
#  CM_obs <- apply(commat_obs[,ausw], 1, function(x) mean(dat$traitmat[x==1]))
#  CM_cor <- apply(z[,ausw], 1, function(x) mean(dat$traitmat[x==1], na.rm = TRUE))
#  CM_true <- apply(dat$z_true, 1, function(x) mean(dat$traitmat[x==1]))
#  
#  # Plot CM along gradient
#  plot(dat$gradient, CM_obs, pch = 16, cex = 0.7, ylim = c(-1.5, 1.5))
#  points(dat$gradient, CM_true, cex = 0.7, col = "orange")
#  points(dat$gradient, CM_cor, pch = "+", cex = 0.7, col = "red")
#  abline(h=0)
#  abline(h=mean(dat$traitmat), lty = 2)
#  
#  
#  ausw <- commat_obs==0 & z == 1
#  plot(dat$gradient, apply(ausw, 1, sum))
#  plot(dat$gradient, apply(ausw, 1, function(x) ))
#  

