# Functions for the book Applied Hierarchical Modeling in Ecology (AHM)
# Marc Kery & Andy Royle, Academic Press, 2016.

#' Function to simulate community occupancy data
#' 
#' Function to simulate community occupancy data with random species effects for
#' psi and p (both including effects of one covariate, 'gradient' for psi and 
#' 'date' for p)
#' 
#' Function simulates data from repeated sampling of a metacommunity according 
#' the model of Dorazio & Royle (JASA, 2005) for type = "det/nondet" (this is 
#' the default) or under the model of Yamaura et al. (2012) for type = "counts".
#' Occupancy probability (psi) or expected abundance (lambda) can be made 
#' dependent on a continuous site covariate 'gradient', while detection 
#' probability can be made dependent an observational covariate 'date'. Both 
#' intercept and slope of the two log- or logistic regressions (for occupancy or
#' expected abundance, respectively, and for detection) are simulated as draws 
#' from a normal distribution with mean and standard deviation that can be 
#' selected using function arguments.
#' 
#' @param nsite Number of sites in the meta-community.
#' @param nrep Number of replicate samples of the meta-community.
#' @param nspec Total number of species in the area that is sampled by these 
#'   sites (i.e. the regional species pool).
#' @param mean.psi Community mean of occupancy probability over all species in 
#'   community (probability scale).
#' @param sig.lpsi Community standard deviation of qlogis(occupancy probability 
#'   intercept).
#' @param mu.FTfilter.lpsi Effect of functional trait on occupancy probability 
#'   (i.e. environmental filtering).
#' @param mu.beta.lpsi Community mean of the effects of 'gradient'covariate on 
#'   occupancy probability.
#' @param sig.beta.lpsi Community standard deviation of the effects of 
#'   'gradient' covariate on occupancy probability.
#' @param mean.p Community mean of detection probability over all species in 
#'   superpopulation (probability scale).
#' @param sig.lp Community standard deviation of qlogis(detection probability 
#'   intercept).
#' @param mu.FTfilter.lp Effect of functional trait on detection probability 
#'   (i.e. detection filtering).
#' @param mu.beta.lp Community mean of the effects of 'date' covariate on 
#'   detection probability.
#' @param sig.beta.lp Community standard deviation of the effects of 
#'   'date'covariate on detection probability.
#'   
#' @return To be written
#'   
#' @examples
#' # Simulate community data
#' set.seed(1234)
#' dat <- simcom(mu.FTfilter.lp = 1.5, mu.FTfilter.lpsi = -0.5)
#' 
#' # Calculate community mean value of trait expression (CM)
#' commat_obs <- apply(dat$y, c(1,2), max)
#' CM_obs <- apply(commat_obs, 1, function(x) mean(dat$traitmat[x==1]))
#' CM_true <- apply(dat$z_true, 1, function(x) mean(dat$traitmat[x==1]))
#' 
#' # Plot CM along gradient
#' plot(dat$gradient, CM_obs, pch = 16, cex = 0.7, ylim = c(-1.5, 1.5))
#' points(dat$gradient, CM_true, cex = 0.7, col = "orange")
#' abline(h=0)
#' abline(h=mean(dat$traitmat), lty = 2)
#' 
#' # Effect of detection filtering along gradient
#' plot(dat$gradient, CM_obs - CM_true, pch = 16, cex = 0.7)
#' abline(h=0)
#' 
#' # Plot detection probability of a species and its trait expression
#' plot(dat$traitmat, dat$P, pch = 16, cex = 0.7, ylim = c(0,1),
#'   xlab = "Expression of functional trait", ylab = "")
#' mtext("Detection probability", 2, line = 3)
#' mtext("(at average gradient)", 2, line = 2.4, cex = 0.8)
#' 
#' @export
#' 
simcom <- function(nsite=50, nrep=3, nspec=100,
                   mean.psi=0.5, sig.lpsi=1, mu.FTfilter.lpsi=0, mu.beta.lpsi=0, sig.beta.lpsi=0,
                   mean.p=0.8, sig.lp=1, mu.FTfilter.lp=0, mu.beta.lp=0, sig.beta.lp=0
                   ) {
  
  # nsite=30; nrep=2; nspec=100;
  # mean.psi=0.5; sig.lpsi=1; mu.FTfilter.lpsi=0; mu.beta.lpsi=0; sig.beta.lpsi=0;
  # mean.p=0.8; sig.lp=1; mu.FTfilter.lp=0; mu.beta.lp=0; sig.beta.lp=0
  # 
  
  # Prepare structures to hold data
  y <- p <- array(NA, c(nsite, nspec, nrep))
  dimnames(y) <- dimnames(p) <-
    list(paste("site", 1:nsite, sep=""), paste("sp", 1:nspec, sep=""), paste("rep", 1:nrep, sep=""))
  z <- psi <- matrix(NA, nsite, nspec)
  dimnames(z) <- dimnames(psi) <- 
    list(paste("site", 1:nsite, sep=""), paste("sp", 1:nspec, sep=""))
  detected.at.all <- rep(NA, nspec)
  
  # Create covariates 'gradient' and 'date'
  gradient <- sort(rnorm(nsite, 0, 1.5))     # Note 'gradient gradient' due to sorting
  date <- matrix(rnorm(nsite * nrep, 0, 1.5), ncol=nrep)
  
  # Create species traits 
  traitmat <- rnorm(nspec, 0, 1.5)
  
  # Draw species-specific intercepts and slopes from their normal distributions
  # Build up linear predictors for occupancy and detection
  a0 <- rnorm(nspec, qlogis(mean.p) + mu.FTfilter.lp * traitmat,  sig.lp)         # detection intercept
  a1 <- rnorm(nspec, mu.beta.lp, sig.beta.lp)               # detection slope on date
  
  # b0 <- qlogis(mean.psi) + mu.FTfilter.lpsi * traitmat      # occupancy intercept
  b0 <- rnorm(nspec, qlogis(mean.psi), sig.lpsi)      # occupancy intercept
  # b1 <- rnorm(nspec, mu.beta.lpsi, sig.beta.lpsi)           # occupancy slope on gradient
  b1 <- rnorm(nspec, mu.beta.lpsi + mu.FTfilter.lpsi * traitmat, sig.beta.lpsi)          # occupancy slope on gradient
  for(k in 1:nspec){
    psi[,k] <- plogis(b0[k] + b1[k] * gradient)
    for(j in 1:nrep){
      p[,k,j] <- plogis(a0[k] + a1[k] * date[,j])
    }
  }
  
  # Distribute species over sites (simulate true state)
  for(k in 1:nspec){
    z[,k] <- rbinom(nsite, 1, psi[,k])
  }
  occurring.in.sample <- apply(z, 2, max) # Presence/absence at study sites
  
  # Measurement of presence/absence (simulate observation)
  for(k in 1:nspec) {
    for(i in 1:nsite){
      for(j in 1:nrep) {
        y[i,k,j] <- rbinom(1, z[i,k], p[i,k,j])
      }
    }
    detected.at.all[k] <- any(y[, k,] > 0)
  }
  
  # Drop species never detected 
  y <- y[,detected.at.all,] 

  # Return data
  list(z_true = z, P = apply(p, 2, mean),
       y = y, detected.at.all = detected.at.all,
       gradient = gradient, date = date,
       traitmat = traitmat)
}

