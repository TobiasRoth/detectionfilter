## ---- eval=FALSE---------------------------------------------------------
#  devtools::install_github("TobiasRoth/detectionfilter")

## ---- echo=TRUE, message=FALSE, warning=FALSE----------------------------
library(detectionfilter)  
library(RColorBrewer)
library(missForest)
library(geometry)
library(mgcv)
library(FD)

## ---- cache=TRUE---------------------------------------------------------
# Prepare data.frame with all data for linear model
d <- traitmat_NA
d$P <- P
d$nobs <- apply(commat, 2, sum)
d$nobs_sd <- as.vector(scale(d$nobs))

# Add for each species the number of occupied plots (and standardize)
d$nocc <- apply(z, 2, sum)
d$nocc_sd <- as.vector(scale(log(d$nocc)))

# Add for each species the average elevation (and standardize)
d$meanel <- apply(z, 2, function(x) mean(plantsBDM$elevation[x==1]))
d$meanel_sd <- as.vector(scale(log(d$meanel+1)))

# Linear model
mod <- lm(qlogis(P) ~ sla + ch + sm + nocc_sd + meanel_sd, data = d)
round(summary(mod)$coef, 3)

## ---- cache=TRUE---------------------------------------------------------
difdays <- plantsBDM$dates$Dat2 - plantsBDM$dates$Dat1
round(summary(lm(difdays ~ I(plantsBDM$elevation/100)))$coef, 3)

## ---- cache=TRUE---------------------------------------------------------

# Species richness SR
SR.obs <- apply(commat, 1, sum)
SR.cor <- apply(z, 1, sum)

# Mean elevation
tmp <- apply(commat, 2, function(x) mean(plantsBDM$elevation[x==1]))
meanel.obs <- apply(commat, 1, function(x) mean(tmp[x==1]))
meanel.cor <- apply(z != commat, 1, function(x) mean(tmp[x==1]))

# Number of occurrences
tmp <- apply(z, 2, sum)
occ.obs <- apply(commat, 1, function(x) mean(tmp[x==1]))
occ.cor <- apply(z != commat, 1, function(x) mean(tmp[x==1]))


## ----fig1, cache=TRUE, fig.width = 10, fig.height = 3.5, echo = FALSE----
# Plot settings
# Plot settings
par(mfrow = c(1, 3), mar = c(2.5,4,3,1.5), oma = c(1,1,1,1))
ele <- 400:2500
xax <- seq(250, 2750, 250)
tcol <- brewer.pal(8, "Set1")
d <- data.frame(SR.obs = SR.obs, SR.cor = SR.cor, prop = SR.obs/SR.cor, meanel.cor=meanel.cor, meanel.obs=meanel.obs,
                occ.obs = occ.obs, occ.cor = occ.cor, elevation = plantsBDM$elevation)

# Proportion of observed species
plot(d$elevation, d$prop, ylim = c(0.7, 1), xlim = c(250,2750),
     pch = 21, cex = 0.8, axes=F, xlab = "", ylab = "", col = "grey50", bg = "grey80")
pred <- predict(gam(prop ~ s(elevation), data = d),
                newdata = data.frame(elevation=ele), type = "response")
points(ele, pred, ty = "l", lty = 1, lwd = 2, col = "black")
axis(side=1, at = xax, labels = rep("", length(xax)), pos = 0.7)
mtext(seq(500,2500,500), 1, at = seq(500,2500,500), cex = 0.7)
axis(side = 2, at = seq(0.7, 1, 0.05), pos = 250, las = 1, labels = rep("", 7))
mtext(seq(0.7, 1, 0.05), 2, at = seq(0.7, 1, 0.05), cex = 0.7, las = 1)
mtext(text = "Community Elevation", side = 1, line = 1, cex = 0.7)
mtext(text = "Proportion of observed species", side = 2, line = 3, cex = 0.7)
mtext(text = "(a) Proportion of observed species\n", side = 3, at = -420, line = 1.5, cex = 0.8, adj = 0)

# Mean elevation
plot(d$elevation, d$meanel.obs, ylim = c(500, 2000), xlim = c(250,2750),
     pch = 21, cex = 0.8, axes=F, xlab = "", ylab = "", col = "grey50", bg = tcol[1])
points(d$elevation, d$meanel.cor, pch = 21, cex = 0.8, col = "grey50", bg = tcol[2])
pred <- predict(gam(meanel.obs ~ s(elevation), data = d),
                newdata = data.frame(elevation=ele), type = "response")
points(ele, pred, ty = "l", lty = 1, lwd = 2, col = tcol[1])
pred <- predict(gam(meanel.cor ~ s(elevation), data = d),
                newdata = data.frame(elevation=ele), type = "response")
points(ele, pred, ty = "l", lty = 1, lwd = 2, col = tcol[2])
axis(side=1, at = xax, labels = rep("", length(xax)), pos = 500)
mtext(seq(500,2500,500), 1, at = seq(500,2500,500), cex = 0.7)
axis(side = 2, at = seq(500, 2000, 250), pos = 250, las = 1, labels = rep("", 7))
mtext(seq(500, 2000, 250), 2, at = seq(500, 2000, 250), cex = 0.7, las = 1)
mtext(text = "Community Elevation", side = 1, line = 1, cex = 0.7)
mtext(text = "Mean elevation of occurrence (m)", side = 2, line = 3, cex = 0.7)
mtext(text = "(b) Community averages of species' mean\nelevation of occurence", side = 3, at = -420, line = 1.5, cex = 0.8, adj = 0)
legend(1500, 730, c("observed communities", "overlooked species only"),
       bty = "n", cex = 0.8, lty = 1, y.intersp=0.8, col = tcol[1:2], pch = 16)

# Number of occurrences
plot(d$elevation, d$occ.obs, ylim = c(50, 250), xlim = c(250,2750),
     pch = 21, cex = 0.8, axes=F, xlab = "", ylab = "", col = "grey50", bg = tcol[1])
points(d$elevation, d$occ.cor, pch = 21, cex = 0.8, col = "grey50", bg = tcol[2])
pred <- predict(gam(occ.obs ~ s(elevation), data = d),
                newdata = data.frame(elevation=ele), type = "response")
points(ele, pred, ty = "l", lty = 1, lwd = 2, col = tcol[1])
pred <- predict(gam(occ.cor ~ s(elevation), data = d),
                newdata = data.frame(elevation=ele), type = "response")
points(ele, pred, ty = "l", lty = 1, lwd = 2, col = tcol[2])
axis(side=1, at = xax, labels = rep("", length(xax)), pos = 50)
mtext(seq(500,2500,500), 1, at = seq(500,2500,500), cex = 0.7)
axis(side = 2, at = seq(50, 250, 50), pos = 250, las = 1, labels = rep("", 5))
mtext(seq(50, 250, 50), 2, at = seq(50, 250, 50), cex = 0.7, las = 1)
mtext(text = "Community Elevation", side = 1, line = 1, cex = 0.7)
mtext(text = "Number of occurences", side = 2, line = 3, cex = 0.7)
mtext(text = "(c) Community averages of species occurrences\n", side = 3, at = -420, line = 1.5, cex = 0.8, adj = 0)
legend(300, 80, c("observed communities", "overlooked species only"),
       bty = "n", cex = 0.8, lty = 1, y.intersp=0.8, col = tcol[1:2], pch = 16)

## ---- cache=TRUE---------------------------------------------------------
# Species that remain most often undetected at communities below 750m
undet <- apply(z[plantsBDM$elevation < 750, ] != commat[plantsBDM$elevation < 750, ], 2, sum)
row.names(traitmat)[order(undet, decreasing = TRUE)][1:2]

# Species that remain most often undetected at communities between 750m and 1250m
undet <- apply(z[plantsBDM$elevation > 750 & plantsBDM$elevation < 1250, ] !=
                 commat[plantsBDM$elevation > 750 & plantsBDM$elevation < 1250, ], 2, sum)
row.names(traitmat)[order(undet, decreasing = TRUE)][1]

## ---- cache=TRUE---------------------------------------------------------
# Specific leaf area (SLA)
f.sla <- function(x) mean(traitmat$sla[as.logical(x)])
CM.sla.obs <- apply(commat, 1, f.sla)
CM.sla.cor <- apply(z, 1, f.sla)
rel <- 100 * abs(lm(CM.sla.cor ~ plantsBDM$elevation)$coef[2])
dif.sla <- abs(CM.sla.obs - CM.sla.cor) > rel

# Canopy height (CH)
f.ch <- function(x) mean(traitmat$ch[as.logical(x)])
CM.ch.obs <- apply(commat, 1, f.ch)
CM.ch.cor <- apply(z, 1, f.ch)
rel <- 100 * abs(lm(CM.ch.cor ~ plantsBDM$elevation)$coef[2])
dif.ch <- abs(CM.ch.obs - CM.ch.cor) > rel

# Seed mass (SM)
f.sm <- function(x) mean(traitmat$sm[as.logical(x)])
CM.sm.obs <- apply(commat, 1, f.sm)
CM.sm.cor <- apply(z, 1, f.sm)
rel <- 100 * abs(lm(CM.sm.cor ~ plantsBDM$elevation)$coef[2])
dif.sm <- abs(CM.sm.obs - CM.sm.cor) > rel

## ------------------------------------------------------------------------
f.plotFD <- function(d, title, tylab, ymin = -0.8, leg.cex = 0.8, lxtext = 3,
                     ymax = 0.6, tticks = 0.2, yleg = c(-0.475, -0.59), 
                     xleg = c(350, 300), confint = FALSE) {
  ele <- 400:2500
  xax <- seq(250, 2750, 250)
  if(!confint) tcol <- brewer.pal(8, "Set1")
  if(confint) {
    tcol <- brewer.pal(4, "Paired")
    tcol <- paste0(tcol, c("50", "99", "50", 99))
  }
  plot(NA, ylim = c(ymin,ymax), xlim = c(250,2750), axes=F, xlab = "", ylab = "")
  farbe <- rep("grey80", nplots)
  farbe[d$dif & d$obs < d$cor] <- tcol[1]
  farbe[d$dif & d$obs > d$cor] <- tcol[2]
  if(!confint) {
    points(d$elevation, d$obs, pch = 21, cex = 0.5, col = "grey50")
    points(d$elevation, d$obs, pch = 16, cex = 0.5, col = farbe)
  }
  if(confint) {
    nsim <- 1000
    tmp <- array(0, dim = c(length(ele), nsim))
    for(s in 1:nsim) {
      tmp[, s] <- predict(gam(obs ~ s(elevation), 
                              data = d[sample(1:nrow(d), nrow(d), replace = TRUE),]), 
                          newdata = data.frame(elevation=ele), type = "response")
    }
    polygon(c(ele, rev(ele)), 
            c(apply(tmp, 1, quantile, probs = 0.025), 
              rev(apply(tmp, 1, quantile, probs = 0.975))), col = tcol[1], border = tcol[1])
    tmp <- array(0, dim = c(length(ele), nsim))
    for(s in 1:nsim) {
      tmp[, s] <- predict(gam(cor ~ s(elevation), 
                              data = d[sample(1:nrow(d), nrow(d), replace = TRUE),]),
                          newdata = data.frame(elevation=ele), type = "response")
    }
    polygon(c(ele, rev(ele)), 
            c(apply(tmp, 1, quantile, probs = 0.025), 
              rev(apply(tmp, 1, quantile, probs = 0.975))), col = tcol[3], border = tcol[3])
  }
  pred <- predict(gam(obs ~ s(elevation), data = d), 
                  newdata = data.frame(elevation=ele), type = "response")
  points(ele, pred, ty = "l", lty = 2, col = ifelse(confint, tcol[2], "black"))
  pred <- predict(gam(cor ~ s(elevation), data = d), 
                  newdata = data.frame(elevation=ele), type = "response")
  points(ele, pred, ty = "l", col = ifelse(confint, tcol[4], "black"))
  axis(side=1, at = xax, labels = rep("", length(xax)), pos=ymin)
  mtext(seq(500,2500,500), 1, at = seq(500,2500,500), cex = 0.7)
  axis(side = 2, at = seq(ymin, ymax, tticks), pos = 250, las = 1, 
       labels = rep("", length(seq(ymin, ymax, tticks))))
  mtext(round(seq(ymin, ymax, tticks), 2), 2, at = seq(ymin, ymax, tticks), 
        cex = 0.7, las = 1)
  mtext(text = "Community Elevation", side = 1, line = 1, cex = 0.7)
  mtext(text = tylab, side = 2, line = lxtext, cex = 0.7)
  mtext(text = title, side = 3, at = -420, line = 1.5, cex = 0.8, adj = 0)
  if(!confint) {
  legend(xleg[1], yleg[1], c("     underestimated values", 
                             "     overestimated values"), col = tcol[1:2], 
         bty = "n", cex = leg.cex, pch = c(16,16), pt.cex = 1, y.intersp=0.8)
  }
  legend(xleg[2], yleg[2], c("prediction from observed communities", 
                             "detection-corrected prediction"), 
         bty = "n", cex = leg.cex, lty = c(2,1), y.intersp=0.8,
         col = ifelse(confint, tcol[c(2,4)], c("black", "black")))
}

## ----fig2.1, cache=TRUE, fig.width = 10, fig.height = 3.5, echo = FALSE----
# Plot settings 
par(mfrow = c(1, 3), mar = c(2.5,4,3,1.5), oma = c(1,1,1,1)) 

# Specific leaf area (SLA)
d <- data.frame(cbind(obs=CM.sla.obs, cor=CM.sla.cor, dif = dif.sla, elevation = plantsBDM$elevation))
f.plotFD(d, "(a) Specific leaf area", "Normalized log-transformed specific leaf area")

# Canopy height (CH)
d <- data.frame(cbind(obs=CM.ch.obs, cor=CM.ch.cor, dif = dif.ch, elevation = plantsBDM$elevation))
f.plotFD(d, "(b) Canopy height", "Normalized log-transformed canopy height")

# Seed mass (SM)
d <- data.frame(cbind(obs=CM.sm.obs, cor=CM.sm.cor, dif = dif.sm, elevation = plantsBDM$elevation))
f.plotFD(d, "(c) Seed mass", "Normalized log-transformed seed mass")

## ----fig2.2, cache=TRUE, fig.width = 10, fig.height = 3.5, echo = FALSE----
# Plot settings
par(mfrow = c(1, 3), mar = c(2.5,4,3,1.5), oma = c(1,1,1,1)) 

# Specific leaf area (SLA)
d <- data.frame(cbind(obs=CM.sla.obs, cor=CM.sla.cor, dif = dif.sla, elevation = plantsBDM$elevation))
f.plotFD(d, "(a) Specific leaf area", "Normalized log-transformed specific leaf area", confint = TRUE)

# Canopy height (CH)
d <- data.frame(cbind(obs=CM.ch.obs, cor=CM.ch.cor, dif = dif.ch, elevation = plantsBDM$elevation))
f.plotFD(d, "(b) Canopy height", "Normalized log-transformed canopy height", confint = TRUE)

# Seed mass (SM)
d <- data.frame(cbind(obs=CM.sm.obs, cor=CM.sm.cor, dif = dif.sm, elevation = plantsBDM$elevation))
f.plotFD(d, "(c) Seed mass", "Normalized log-transformed seed mass", confint = TRUE)

## ---- cache=TRUE---------------------------------------------------------
# Functional trait space (FRic)
f.FRic <- function(x) convhulln(traitmat[as.logical(x),], "FA")$vol
FRic.obs <- apply(commat, 1, f.FRic)
FRic.cor <- apply(z, 1, f.FRic)
rel <- 100 * abs(lm(FRic.obs ~ plantsBDM$elevation)$coef[2])
FRic.dif <- abs(FRic.obs - FRic.cor) > rel

## Calculate distance matrix from trait matrix 
dist <- as.matrix(dist(traitmat))

# Mean nearest neighbout distance (mnnd)
f.mnnd <- function(x) {
  sample.dis <- dist[as.logical(x),as.logical(x)]
  mean(apply(sample.dis, 1, function(x) min(x[x>0])))
}
mnnd.obs <- apply(commat, 1, f.mnnd)
mnnd.cor <- apply(z, 1, f.mnnd)
rel <- 100 * abs(lm(mnnd.obs ~ plantsBDM$elevation)$coef[2])
mnnd.dif <- abs(mnnd.obs - mnnd.cor) > rel

# Species richness SR
SR.obs <- apply(commat, 1, sum)
SR.cor <- apply(z, 1, sum)
rel <- 100 * abs(lm(SR.obs ~ plantsBDM$elevation)$coef[2])
SR.dif <- abs(SR.obs - SR.cor) > rel

## ----fig3.1, cache=TRUE, fig.width = 10, fig.height = 3.5, echo = FALSE----
# Plot settings
par(mfrow = c(1, 3), mar = c(2.5,4,3,1.5), oma = c(1,1,1,1))

# Functional richness
d <- data.frame(cbind(obs=FRic.obs, cor=FRic.cor, dif = FRic.dif, elevation = plantsBDM$elevation))
f.plotFD(d, "(a) Functional richness", "Convex hull volume (FRic)", 
         ymin = 0, ymax = 200, tticks = 25, yleg = c(40, 22))

# Functional packing
d <- data.frame(cbind(obs = mnnd.obs, cor = mnnd.cor, dif = mnnd.dif, elevation = plantsBDM$elevation))
f.plotFD(d, "(b) Functional packing", "Mean nearest neighbour distance", 
         ymin = 0.2, ymax = 0.5, tticks = 0.05, yleg = c(0.26, 0.235))

# Species richness
d <- data.frame(cbind(obs=SR.obs, cor=SR.cor, dif = SR.dif, elevation = plantsBDM$elevation))
f.plotFD(d, "(c) Taxonomic richness", "Number of species", 
         ymin = 0, ymax = 500, tticks = 100, yleg = c(100, 55))

## ----fig3.2, cache=TRUE, fig.width = 10, fig.height = 3.5, echo = FALSE----
# Plot settings 
par(mfrow = c(1, 3), mar = c(2.5,4,3,1.5), oma = c(1,1,1,1))

# Functional richness
d <- data.frame(cbind(obs=FRic.obs, cor=FRic.cor, dif = FRic.dif, elevation = plantsBDM$elevation))
f.plotFD(d, "(a) Functional richness", "Convex hull volume (FRic)", 
         ymin = 0, ymax = 200, tticks = 25, yleg = c(40, 22), confint = TRUE)

# Functional packing
d <- data.frame(cbind(obs = mnnd.obs, cor = mnnd.cor, dif = mnnd.dif, elevation = plantsBDM$elevation))
f.plotFD(d, "(b) Functional packing", "Mean nearest neighbour distance", 
         ymin = 0.2, ymax = 0.5, tticks = 0.05, yleg = c(0.26, 0.235), confint = TRUE)

# Species richness
d <- data.frame(cbind(obs=SR.obs, cor=SR.cor, dif = SR.dif, elevation = plantsBDM$elevation))
f.plotFD(d, "(c) Taxonomic richness", "Number of species", 
         ymin = 0, ymax = 500, tticks = 100, yleg = c(100, 55), confint = TRUE)

## ---- cache=TRUE---------------------------------------------------------
nsim <- 10

# Function to get a random sample of 
tsam <- function(x) {
  z[x, sample(which(z[x,] == 1), SR.cor[x] - SR.obs[x], replace = FALSE)] <- 0
  z[x,]
}

# Run simulation and calculate mean over all simulations
simFRic <- simmnnd <- array(0, dim = c(nrow(commat), nsim))
for(s in 1:nsim) {
  new.z <- t(sapply(1:nrow(commat), tsam))
  simFRic[, s] <- apply(new.z, 1, f.FRic)
  simmnnd[, s] <- apply(new.z, 1, f.mnnd)
}

FRic.cor <- apply(simFRic, 1, mean)
mnnd.cor <- apply(simmnnd, 1, mean)

## ----fig3.3, cache = TRUE, fig.width = 10, fig.height = 3.5, echo = FALSE----
# Plot settings 
par(mfrow = c(1, 3), mar = c(2.5,4,3,1.5), oma = c(1,1,1,1))

# Functional richness
d <- data.frame(cbind(obs=FRic.obs, cor=FRic.cor, dif = FRic.dif, elevation = plantsBDM$elevation))
f.plotFD(d, "(a) Functional richness", "Convex hull volume (FRic)", 
         ymin = 0, ymax = 200, tticks = 25, yleg = c(40, 22), confint = TRUE)

# Functional packing
d <- data.frame(cbind(obs = mnnd.obs, cor = mnnd.cor, dif = mnnd.dif, elevation = plantsBDM$elevation))
f.plotFD(d, "(b) Functional packing", "Mean nearest neighbour distance", 
         ymin = 0.2, ymax = 0.5, tticks = 0.05, yleg = c(0.26, 0.235), confint = TRUE)

## ------------------------------------------------------------------------
# Specific leaf area (SLA)
f <- function(x) mean(traitmat$sla[as.logical(x)])
used <- apply(commat, 1, f)
f <- function(x) mean(traitmat_all$sla[as.logical(x)])
allobs <- apply(commat_obs, 1, f)
cor(used, allobs)

# Canopy height (CH)
f <- function(x) mean(traitmat$ch[as.logical(x)])
used <- apply(commat, 1, f)
f <- function(x) mean(traitmat_all$ch[as.logical(x)])
allobs <- apply(commat_obs, 1, f)
cor(used, allobs)

# Seed mass (SM)
f <- function(x) mean(traitmat$sm[as.logical(x)])
used <- apply(commat, 1, f)
f <- function(x) mean(traitmat_all$sm[as.logical(x)])
allobs <- apply(commat_obs, 1, f)
cor(used, allobs)

# Functional trait space (FRic)
f <- function(x) convhulln(traitmat[as.logical(x),], "FA")$vol
used <- apply(commat, 1, f)
f <- function(x) convhulln(traitmat_all[as.logical(x),], "FA")$vol
allobs <- apply(commat_obs, 1, f)
cor(used, allobs)

## Calculate distance matrix from complete trait matrix 
dist_all <- as.matrix(dist(traitmat_all))

# Mean nearest neighbout distance (mnnd)
f <- function(x) {
  sample.dis <- dist[as.logical(x),as.logical(x)]
  mean(apply(sample.dis, 1, function(x) min(x[x>0])))
}
used <- apply(commat, 1, f)
f <- function(x) {
  sample.dis <- dist_all[as.logical(x),as.logical(x)]
  mean(apply(sample.dis, 1, function(x) min(x[x>0])))
}
allobs <- apply(commat_obs, 1, f)
cor(used, allobs)

# Species richness SR
used <- apply(commat, 1, sum)
allobs <- apply(commat_obs, 1, sum)
cor(used, allobs)

## ---- cache=TRUE---------------------------------------------------------
# Specific leaf area (SLA)
f.sla <- function(x) mean(traitmat$sla[as.logical(x)])
CM.sla.obs <- apply(commat, 1, f.sla)
f.sla <- function(x) mean(traitmat_all$sla[as.logical(x)])
CM.sla.cor <- apply(commat_obs, 1, f.sla)
rel <- 100 * abs(lm(CM.sla.cor ~ plantsBDM$elevation)$coef[2])
dif.sla <- abs(CM.sla.obs - CM.sla.cor) > rel

# Canopy height (CH)
f.ch <- function(x) mean(traitmat$ch[as.logical(x)])
CM.ch.obs <- apply(commat, 1, f.ch)
f.ch <- function(x) mean(traitmat_all$ch[as.logical(x)])
CM.ch.cor <- apply(commat_obs, 1, f.ch)
rel <- 100 * abs(lm(CM.ch.cor ~ plantsBDM$elevation)$coef[2])
dif.ch <- abs(CM.ch.obs - CM.ch.cor) > rel

# Seed mass (SM)
f.sm <- function(x) mean(traitmat$sm[as.logical(x)])
CM.sm.obs <- apply(commat, 1, f.sm)
f.sm <- function(x) mean(traitmat_all$sm[as.logical(x)])
CM.sm.cor <- apply(commat_obs, 1, f.sm)
rel <- 100 * abs(lm(CM.sm.cor ~ plantsBDM$elevation)$coef[2])
dif.sm <- abs(CM.sm.obs - CM.sm.cor) > rel

## ---- cache=TRUE---------------------------------------------------------
# Functional trait space (FRic)
f.FRic <- function(x) convhulln(traitmat[as.logical(x),], "FA")$vol
FRic.obs <- apply(commat, 1, f.FRic)
f.FRic <- function(x) convhulln(traitmat_all[as.logical(x),], "FA")$vol
FRic.cor <- apply(commat_obs, 1, f.FRic)
rel <- 100 * abs(lm(FRic.obs ~ plantsBDM$elevation)$coef[2])
FRic.dif <- abs(FRic.obs - FRic.cor) > rel

## Calculate distance matrix from trait matrix 
dist <- as.matrix(dist(traitmat))

# Mean nearest neighbout distance (mnnd)
dist <- as.matrix(dist(traitmat))
f.mnnd <- function(x) {
  sample.dis <- dist[as.logical(x),as.logical(x)]
  mean(apply(sample.dis, 1, function(x) min(x[x>0])))
}
mnnd.obs <- apply(commat, 1, f.mnnd)
dist <- as.matrix(dist(traitmat_all))
f.mnnd <- function(x) {
  sample.dis <- dist[as.logical(x),as.logical(x)]
  mean(apply(sample.dis, 1, function(x) min(x[x>0])))
}
mnnd.cor <- apply(commat_obs, 1, f.mnnd)
rel <- 100 * abs(lm(mnnd.obs ~ plantsBDM$elevation)$coef[2])
mnnd.dif <- abs(mnnd.obs - mnnd.cor) > rel

# Species richness SR
SR.obs <- apply(commat, 1, sum)
SR.cor <- apply(z, 1, sum)
rel <- 100 * abs(lm(SR.obs ~ plantsBDM$elevation)$coef[2])
SR.dif <- abs(SR.obs - SR.cor) > rel

