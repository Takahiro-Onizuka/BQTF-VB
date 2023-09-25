# location function
# Faulkner github



# ----- Piecewise Constant -----


piecewiseConst_traj <- function(x,fv, bp){
  
  num1 <- length(x)
  num2 <- length(bp)
  out <- rep(NA,num1)
      
  for (j in 1:num2) {
    a <- bp[j]
    b <- bp[j+1]
      
    for (i in 1:num1) {
      if(a < i & i <= b){
        out[i] <- fv[j]
      }
    }
  }
  out
}











# ----- Varying smoothness (Boom-bust trajectory) -----

boomBust_traj <- function(x, mf, sf,soff, sscale, snht, bloc, bscale, bht, tmx=NULL){
  if (!is.null(tmx)) {
    if (length(x)==1) {
      if(x > tmx) x <- tmx
    }
    if (length(x)>1) x[x>=tmx] <- tmx
  }	 
  yg <- snht*sin((soff-x)/sscale) + bht*exp(-((x-bloc)^2)/bscale)
  out <- mf + sf*yg
  out
}








# ----- Smooth trend GP -----

expCovih <- function(ssv, pvec){
  # plist is list of parms {sigmaf2, ll}
  sigmaf2 <- pvec[1]
  ll <- pvec[2]
  nn <- length(ssv)
  cmat <- matrix(0, nn, nn)
  for (i in 1:nn){
    for (j in 1:nn){
      cmat[i,j] <- sigmaf2*exp(-0.5*((i-j)/ll)^2)
    }
  }
  cmat
}



Smooth_GP <- function(x, mu, sigma, rho){
  
  pvec <- rep(NA,2)
  pvec[1] <- sigma
  pvec[2] <- rho
  
  Sig <- expCovih(x,pvec)
  
  library(MASS)
  # Generate a number of functions from the process
  values <- mvrnorm(1, rep(mu, length(x)), Sig)
  
  return(values)
}

#set.seed(29)
#th0 <- Smooth_GP((1:100)/100, 2, 1, 10)
#plot(th0, type="l")
