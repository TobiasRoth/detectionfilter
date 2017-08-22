## ---- eval=FALSE---------------------------------------------------------
#  devtools::install_github("TobiasRoth/detectionfilter")

## ---- echo=TRUE, message=FALSE, warning=FALSE----------------------------
library(detectionfilter)  
library(RColorBrewer)
library(missForest)
library(geometry)
library(mgcv)
library(FD)

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

