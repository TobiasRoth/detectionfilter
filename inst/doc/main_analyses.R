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
f.plotFD <- function(d, title, tylab, ymin = -0.75, leg.cex = 0.8, lxtext = 3,
                     ymax = 0.5, tticks = 0.25, yleg = c(-0.475, -0.59), 
                     xleg = c(350, 300)) {
  ele <- 400:2500
  xax <- seq(250, 2750, 250)
  tcol <- brewer.pal(8, "Set1")
  plot(NA, ylim = c(ymin,ymax), xlim = c(250,2750), axes=F, xlab = "", ylab = "")
  farbe <- rep("grey80", nplots)
  farbe[d$dif & d$obs < d$cor] <- tcol[1]
  farbe[d$dif & d$obs > d$cor] <- tcol[2]
  points(d$elevation, d$obs, pch = 21, cex = 0.5, col = "grey50")
  points(d$elevation, d$obs, pch = 16, cex = 0.5, col = farbe)
  pred <- predict(gam(obs ~ s(elevation), data = d), 
                  newdata = data.frame(elevation=ele), type = "response")
  points(ele, pred, ty = "l", lty = 2)
  pred <- predict(gam(cor ~ s(elevation), data = d), 
                  newdata = data.frame(elevation=ele), type = "response")
  points(ele, pred, ty = "l")
  axis(side=1, at = xax, labels = rep("", length(xax)), pos=ymin)
  mtext(seq(500,2500,500), 1, at = seq(500,2500,500), cex = 0.7)
  axis(side = 2, at = seq(ymin, ymax, tticks), pos = 250, las = 1, 
       labels = rep("", length(seq(ymin, ymax, tticks))))
  mtext(round(seq(ymin, ymax, tticks), 2), 2, at = seq(ymin, ymax, tticks), 
        cex = 0.7, las = 1)
  mtext(text = "Community Elevation", side = 1, line = 1, cex = 0.7)
  mtext(text = tylab, side = 2, line = lxtext, cex = 0.7)
  mtext(text = title, side = 3, at = -420, line = 1.5, cex = 0.8, adj = 0)
  legend(xleg[1], yleg[1], c("     underestimated values", 
                             "     overestimated values"), col = tcol[1:2], 
         bty = "n", cex = leg.cex, pch = c(16,16), pt.cex = 1, y.intersp=0.8)
  legend(xleg[2], yleg[2], c("prediction from observed communities", 
                             "detection corrected prediction"), 
         bty = "n", cex = leg.cex, lty = c(2,1), y.intersp=0.8)
}

