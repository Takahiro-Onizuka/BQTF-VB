source("BQTF-function.R")
source("true-function.R")

# function -----------
quant=function(x){ quantile(x, prob=c(0.025,0.975)) }
# --------------------

## setting
n <- 100
pp <- 0.5  # quantile level
mc <- 7500 # MCMC sample size
th <- 10
thin <- seq(from=th, to=mc, by=th) # thinning 

## data generating process
TrueFun <- "VS"
set.seed(1)
X <- 1:n
if(TrueFun=="PC"){
  th0 <- piecewiseConst_traj(1:n, fv=c(2.5, 1,3.5,1.5), bp=c(0,20,40,60,100))
  TForder <- 0
}else if(TrueFun=="VS"){
  th0 <- boomBust_traj(4*(1:n)/n, mf=2,sf=1 ,soff=2, sscale=-1, snht=1, bloc=2, bscale=1/30, bht=2, tmx=NULL)
  TForder <- 1
}
Y <- th0 + rnorm(n,0,(1+(X/n)^2)/4)

# true quantile trend
True <- th0 + qnorm(pp, 0, (1+((1:n)/n)^2)/4)

# -----------------#
#     fitting      #
# -----------------#

# Fitting of MCMC-HS
fit <- BQTF.MCMC.HS.multi(X=X, Y=Y, TForder=TForder, mc=mc, pp=pp)
Est <- apply(fit[thin,], 2, mean)
CI <- apply(fit[thin,], 2, quant)

plot(Y, col=1, ylim=c(0,5), main="MCMC-HS")
points(True, type="l", col=1, ylim=c(-2,2))
points(th0, type="l", col=1, ylim=c(-2,2),lty=2, lwd=2)
points(Est, type="l", col=2, ylim=c(-2,2))
polygon(c(X,rev(X)), c(CI[1,],rev(CI[2,])), col=adjustcolor(2, alpha=0.4), border=NA)





# Fitting of VB-HS
fit <- BQTF.VB.HS.multi(X=X, Y=Y, TForder=TForder, pp=pp)
Est <- fit$sm
CI <- matrix(NA, 2, n)
CI[1,] <- qnorm(0.025, fit$sm, sqrt(fit$var))
CI[2,] <- qnorm(0.975, fit$sm, sqrt(fit$var))

plot(Y, col=1, ylim=c(0,5), main="VB-HS")
points(True, type="l", col=1, ylim=c(-2,2))
points(th0, type="l", col=1, ylim=c(-2,2),lty=2, lwd=2)
points(Est, type="l", col=2, ylim=c(-2,2))
polygon(c(X,rev(X)), c(CI[1,],rev(CI[2,])), col=adjustcolor(2, alpha=0.4), border=NA)





# Fitting of CVB-HS
fit <- BQTF.CVB.multi(X=X, Y=Y, TForder=TForder, pp=pp, prior="HS", B.size=100)
Est <- fit$sm
CI <- matrix(NA, 2, n)
CI[1,] <- qnorm(0.025, fit$sm, sqrt(fit$var))
CI[2,] <- qnorm(0.975, fit$sm, sqrt(fit$var))

plot(Y, col=1, ylim=c(0,5), main="CVB-HS")
points(True, type="l", col=1, ylim=c(-2,2))
points(th0, type="l", col=1, ylim=c(-2,2),lty=2, lwd=2)
points(Est, type="l", col=2, ylim=c(-2,2))
polygon(c(X,rev(X)), c(CI[1,],rev(CI[2,])), col=adjustcolor(2, alpha=0.4), border=NA)



