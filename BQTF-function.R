library(MASS)
library(GIGrvg)
library(invgamma)
library(Bessel)
library(MCMCpack)




##---------------------------------------------##
##       BQTF-MCMC (multiple observation)      ##
##---------------------------------------------##




BQTF.MCMC.HS.multi <- function(X, Y, TForder, pp=0.5, mc=1000, sparse=F, sig.init=c(1,1)){
  
  Data <- data.frame(cbind(X,Y))
  colnames(Data) <- c("X","Y")
  Data <- Data[order(Data$X),]
  X <- Data$X
  Y <- Data$Y
  location <- unique(Data$X)
  n <- length(location)
  N <- length(Y)
  Ni <- rep(NA, n)
  y <- rep(NA, n)
  ymax <- rep(NA, n)
  num <- rep(0, n+1)
  for (i in 1:n) {
    data <- subset(Data, Data$X == location[i])
    Ni[i] <- length(data$Y)
    y[i] <- sum(data$Y)
    num[i+1] <- num[i]+Ni[i]
    ymax[i] <- max(data$Y)
  }
  
  ## settings
  #n <- length(y)  
  psi <- (1-2*pp)/(pp*(1-pp))
  t2 <- 2/(pp*(1-pp))
  
  # difference operator
  D1 <- matrix(0, n-1,n)
  for (i in 1:n-1) {
    for (j in 1:n) {
      if(i == j){
        D1[i,j] <- 1
      }
      if(i+1 == j){
        D1[i,j] <- -1
      }
    }
  }
  if(TForder == 0){
    DD <- D1
    D <- rbind(c(1,rep(0, n-1)),DD)
  }
  if(TForder == 1){
    DD <- D1[1:(n-2),1:(n-1)]%*%diag(TForder/diff(location))%*%D1
    D <- rbind(t(matrix(c(1,rep(0, n-1), 
                          0,1,rep(0, n-2)), n, 2)),DD)
  }
  if(TForder == 2){
    D2 <- D1[1:(n-2),1:(n-1)]%*%diag(TForder/diff(location))%*%D1
    DD <- D1[1:(n-3),1:(n-2)]%*%diag(TForder/diff(location,lag=2))%*%D2
    D <- rbind(t(matrix(c(1,rep(0, n-1),
                          0,1,rep(0, n-2),
                          0,0, 1, rep(0, n-3)), n, 3)),DD)
  }
  
  TD <- t(D)
  if(sparse){
    TD <- as(TD, "sparseMatrix")
    D <- as(D, "sparseMatrix")
  }
  
  
  ## MCMC box
  th.pos <- matrix(NA, mc, n)
  sig.pos <- c()
  
  ## initial value
  th <- y
  sig2 <- 1
  IGa_w <- rep(1, TForder+1) 
  IGb_w <- rep(1, TForder+1) 
  ww <- rep(1, n)   # local parameter in HS
  lower <- 10^(-6)   # lower bound
  upper <- 10^6   # upper bound 
  tau2 <- 1   # global parameter in HS
  tau2vec <- rep(1,n)
  nu <- rep(1, n-TForder-1)   # latent parameter in HS
  xi <- 1   # latent parameter in HS
  vv <- rep(1, N)   # latent parameter in asymmetric Laplace,Exp(sig2)
  IGa <- sig.init[1]   # parameter of sig2 prior
  IGb <- sig.init[2]   # parameter of sig2 prior
  
  ## MCMC replications
  for(k in 1:mc){
    
    # theta
    cc <- rep(0, n)
    cc <- as.vector(sapply(split(vv^(-1), X), sum))/t2
    A <- solve(TD%*%diag(1/(ww*tau2vec))%*%D + diag(cc) )
    A <- (t(A) + A)/2
    B <- as.vector(sapply(split((Y-psi*vv)/vv, X), sum))/t2
    th <- mvrnorm(1, A%*%B, sig2*A)
    th.pos[k,] <- th
    
    ## asymmetric Laplace 
    # sig2
    eta <- as.vector(D%*%th)
    bb <- rep(0, N)
    for (i in 2:(n+1)) {
      bb[(num[i-1]+1):num[i]] <- Y[(num[i-1]+1):num[i]]-th[i-1]
    }
    cc <- sum(((bb-psi*vv)^2)/vv)/(2*t2) + sum(eta^2/(ww*tau2vec))/2 + sum(vv) + IGb
    sig2 <- rinvgamma(1, (n+3*N)/2+IGa, cc)   
    sig.pos[k] <- sig2
    
    # vv
    bb2 <- bb^2/(t2*sig2)
    cc <- 2/sig2 + psi^2/(t2*sig2)
    vv <- tapply(bb2, 1:N, rgig, n=1, lambda=1/2, psi=cc)
    vv[vv<lower] <- lower
    
    ## HS latent variables
    # ww
    #eta <- as.vector(D%*%th)
    bb <- eta^2/(2*tau2vec*sig2)
    ww[1:(TForder+1)] <- rinvgamma(TForder+1, 1/2+IGa_w, bb[1:(TForder+1)]+IGb_w)
    ww[(TForder+2):n] <- rinvgamma(n-TForder-1, 1, bb[(TForder+2):n]+1/nu)
    ww[ww<lower] <- lower
    ww[ww>upper] <- upper
    
    # nu
    nu <- rinvgamma(n-TForder-1, 1/2, 1+1/ww[(TForder+2):n])
    
    # tau2
    bbb <- sum(eta[(TForder+2):n]^2/(ww[(TForder+2):n]*sig2))/2 + 1/xi
    tau2 <- min(max(rinvgamma(1, (n-TForder)/2, bbb), 10^(-5)), 10^5)   
    tau2vec <- c(rep(1,TForder+1), rep(tau2, n-TForder-1))
    
    
    # xi
    xi <- rinvgamma(1, 1/2, 1+1/tau2)
    
    # print 
    #if(round(k/200)==(k/200)){ print(k) }
  }
  
  # output
  return(th.pos)
}




BQTF.MCMC.Lap.multi <- function(X, Y, TForder, pp=0.5, mc=1000, sparse=F, sig.init=c(1,1)){
  
  Data <- data.frame(cbind(X,Y))
  colnames(Data) <- c("X","Y")
  Data <- Data[order(Data$X),]
  X <- Data$X
  Y <- Data$Y
  location <- unique(Data$X)
  n <- length(location)
  N <- length(Y)
  Ni <- rep(NA, n)
  y <- rep(NA, n)
  ymax <- rep(NA, n)
  num <- rep(0, n+1)
  for (i in 1:n) {
    data <- subset(Data, Data$X == location[i])
    Ni[i] <- length(data$Y)
    y[i] <- sum(data$Y)
    num[i+1] <- num[i]+Ni[i]
    ymax[i] <- max(data$Y)
  }
  
  ## settings
  #n <- length(y)  
  psi <- (1-2*pp)/(pp*(1-pp))
  t2 <- 2/(pp*(1-pp))
  
  # difference operator
  D1 <- matrix(0, n-1,n)
  for (i in 1:n-1) {
    for (j in 1:n) {
      if(i == j){
        D1[i,j] <- 1
      }
      if(i+1 == j){
        D1[i,j] <- -1
      }
    }
  }
  if(TForder == 0){
    DD <- D1
    D <- rbind(c(1,rep(0, n-1)),DD)
  }
  if(TForder == 1){
    DD <- D1[1:(n-2),1:(n-1)]%*%diag(TForder/diff(location))%*%D1
    D <- rbind(t(matrix(c(1,rep(0, n-1), 
                          0,1,rep(0, n-2)), n, 2)),DD)
  }
  if(TForder == 2){
    D2 <- D1[1:(n-2),1:(n-1)]%*%diag(TForder/diff(location))%*%D1
    DD <- D1[1:(n-3),1:(n-2)]%*%diag(TForder/diff(location,lag=2))%*%D2
    D <- rbind(t(matrix(c(1,rep(0, n-1),
                          0,1,rep(0, n-2),
                          0,0, 1, rep(0, n-3)), n, 3)),DD)
  }
  
  TD <- t(D)
  if(sparse){
    TD <- as(TD, "sparseMatrix")
    D <- as(D, "sparseMatrix")
  }
  
  
  ## MCMC box
  th.pos <- matrix(NA, mc, n)
  sig.pos <- c()
  
  ## initial value
  th <- y
  sig2 <- 1
  IGa_w <- rep(1, TForder+1) 
  IGb_w <- rep(1, TForder+1) 
  ww <- rep(1, n)   # local parameter in Lap
  gam <- 1
  nu <- 1
  lower <- 10^(-6)   # lower bound
  upper <- 10^6   # upper bound 
  vv <- rep(1, N)   # latent parameter in asymmetric Laplace,Exp(sig2)
  IGa <- sig.init[1]   # parameter of sig2 prior
  IGb <- sig.init[2]   # parameter of sig2 prior
  
  ## MCMC replications
  for(k in 1:mc){
    
    # theta
    cc <- rep(0, n)
    cc <- as.vector(sapply(split(vv^(-1), X), sum))/t2
    A <- solve(TD%*%diag(1/(ww))%*%D + diag(cc) )
    A <- (t(A) + A)/2
    B <- as.vector(sapply(split((Y-psi*vv)/vv, X), sum))/t2
    th <- mvrnorm(1, A%*%B, sig2*A)
    th.pos[k,] <- th
    
    ## asymmetric Laplace 
    # sig2
    eta <- as.vector(D%*%th)
    bb <- rep(0, N)
    for (i in 2:(n+1)) {
      bb[(num[i-1]+1):num[i]] <- Y[(num[i-1]+1):num[i]]-th[i-1]
    }
    cc <- sum(((bb-psi*vv)^2)/vv)/(2*t2) + sum(eta^2/(ww))/2 + sum(vv) + IGb
    sig2 <- rinvgamma(1, (n+3*N)/2+IGa, cc)   
    sig.pos[k] <- sig2
    
    # vv
    bb2 <- bb^2/(t2*sig2)
    cc <- 2/sig2 + psi^2/(t2*sig2)
    vv <- tapply(bb2, 1:N, rgig, n=1, lambda=1/2, psi=cc)
    vv[vv<lower] <- lower
    
    # ww
    #eta <- as.vector(D%*%th)
    bb <- eta^2/sig2
    ww[1:(TForder+1)] <- rinvgamma(TForder+1, 1/2+IGa_w, bb[1:(TForder+1)]/2+IGb_w)
    ww[(TForder+2):n] <- tapply(bb[(TForder+2):n], 1:(n-TForder-1), rgig, n=1, lambda=1/2, psi=gam)
    ww[ww<lower] <- lower
    ww[ww>upper] <- upper
    
    # gam
    gam <- rgig(1, lambda=n-TForder-3/2, chi=2/nu, psi=sum(ww[(TForder+2):n]))
    
    # nu
    nu <- rinvgamma(1, 1/2, 1+1/gam)
    
    
    # print 
    #if(round(k/200)==(k/200)){ print(k) }
  }
  
  # output
  return(th.pos)
}
















##------------------------------------------##
##      BQTF-VB (multiple observation)      ##
##------------------------------------------##

BQTF.VB.HS.multi <- function(X, Y, TForder=1, pp=0.5, epsilon=10^(-3), sig.init=c(0.1,0.1)){
  
  Data <- data.frame(cbind(X,Y))
  colnames(Data) <- c("X","Y")
  Data <- Data[order(Data$X),]
  X <- Data$X
  Y <- Data$Y
  location <- unique(Data$X)
  n <- length(location)
  N <- length(Y)
  Ni <- rep(NA, n)
  y <- rep(NA, n)
  ymax <- rep(NA, n)
  num <- rep(0, n+1)
  for (i in 1:n) {
    data <- subset(Data, Data$X == location[i])
    Ni[i] <- length(data$Y)
    y[i] <- sum(data$Y)
    num[i+1] <- num[i]+Ni[i]
    ymax[i] <- max(data$Y)
  }
  
  
  ## settings
  psi <- (1-2*pp)/(pp*(1-pp))
  t2 <- 2/(pp*(1-pp))
  
  # difference operator
  D1 <- matrix(0, n-1,n)
  for (i in 1:n-1) {
    for (j in 1:n) {
      if(i == j){
        D1[i,j] <- 1
      }
      if(i+1 == j){
        D1[i,j] <- -1
      }
    }
  }
  if(TForder == 0){
    DD <- D1
    D <- rbind(c(1,rep(0, n-1)),DD)
  }
  if(TForder == 1){
    DD <- D1[1:(n-2),1:(n-1)]%*%diag(TForder/diff(location))%*%D1
    D <- rbind(t(matrix(c(1,rep(0, n-1), 
                          0,1,rep(0, n-2)), n, 2)),DD)
  }
  if(TForder == 2){
    D2 <- D1[1:(n-2),1:(n-1)]%*%diag(TForder/diff(location))%*%D1
    DD <- D1[1:(n-3),1:(n-2)]%*%diag(TForder/diff(location,lag=2))%*%D2
    D <- rbind(t(matrix(c(1,rep(0, n-1),
                          0,1,rep(0, n-2),
                          0,0, 1, rep(0, n-3)), n, 3)),DD)
  }
  
  TD <- t(D)
  invD <- solve(D)
  
  # setting
  oldtheta <- 1
  oldvar <- diag(rep(1,n))
  rr <- 1
  rrvar <- 1
  oldXX <- XX <- rep(NA, n)
  
  ## box
  th.pos <- rep(NA, n)
  var.th <- rep(NA, n)
  
  
  # initial value
  IGa <- sig.init[1]   # paramater of sig2 prior
  IGb <- sig.init[2]   # paramater of sig2 prior
  d2 <- 1  # parameter of nu prior
  lower <- 10^(-6)   # lower bound
  upper <- 10^6   # upper bound 
  zlower <- 10^(-6)   # lower bound
  zupper <- 10^6   # upper bound
  a_w <- rep(1, TForder+1)
  b_w <- rep(1, TForder+1)
  
  E_theta <- rep(2,n)
  E_thetai2 <- rep(4,n)
  E_invsig2 <- 1
  E_sig2 <- 1
  a_sig2 <- 1
  E_z <- rep(1,N)
  E_invz <- rep(1,N)
  a_z <- rep(1,N)
  b_z <- 1
  E_invw <- rep(1,n)  
  E_invtau2 <- 1  
  invtau2vec <- rep(1,n)
  E_invnu <- rep(1,n-TForder-1)
  E_invxi <- 0.1
  eta2 <- rep(1,n)
  invB <- diag(rep(1,n))
  gigpara <- 0.1
  sum.invz <- rep(1, n)
  sum.yinvz <- as.vector(sapply(split(Y*E_invz, X), sum))
  
  E_invwinvtau2vec <- E_invw*invtau2vec
  
  for(j in 1:100000){
    
    oldXX <- diag(oldvar)
    XX <- diag(invB/E_invsig2)
    
    rr <- max(abs(oldtheta-E_theta))
    rrvar <- max(abs(oldXX-XX))
    
    if(rr < epsilon & rrvar < epsilon){
      th.pos <- E_theta
      var.th <- diag(invB/E_invsig2)
      break
    }else{
      
      oldtheta <- E_theta
      oldvar <- invB/E_invsig2
      
      # theta
      sum.invz <- as.vector(sapply(split(E_invz, X), sum))/t2
      invB <- solve( diag(sum.invz) + TD%*%diag(E_invwinvtau2vec)%*%D )
      invB <- (t(invB) + invB)/2
      C <- as.vector(sapply(split(Y*E_invz-psi, X), sum))/t2
      AA <- invB/E_invsig2 + invB%*%C%*%t(C)%*%invB
      E_theta <- invB%*%C
      E_thetai2 <- diag(AA)
      eta2 <- diag(D%*%AA%*%TD)
      
      # sigma^2
      aa <- bb <- rep(0, N)
      for (i in 2:(n+1)) {
        aa[(num[i-1]+1):num[i]] <- Y[(num[i-1]+1):num[i]]^2-2*Y[(num[i-1]+1):num[i]]*E_theta[i-1]+E_thetai2[i-1]
        bb[(num[i-1]+1):num[i]] <- Y[(num[i-1]+1):num[i]]-E_theta[i-1]
      }
      a_sig2 <- sum(aa*E_invz - 2*psi*bb + psi^2*E_z)/(2*t2) + sum(eta2*E_invw*invtau2vec)/2 + sum(E_z) + IGb
      E_invsig2 <- ((3*N+n)/2+IGa)/a_sig2
      E_sig2 <- a_sig2/((3*N+n)/2+IGa-1)
      
      # z
      a_z <- aa*E_invsig2/t2
      b_z <- (psi^2*E_invsig2/t2+2*E_invsig2)
      sqab <- sqrt(a_z*b_z)
      E_z[sqab > gigpara] <- as.vector(sqrt(a_z[sqab > gigpara])*BesselK(sqab[sqab > gigpara], 3/2, expon.scaled = FALSE)/
                                         (sqrt(b_z)*BesselK(sqab[sqab > gigpara], 1/2, expon.scaled = FALSE)))
      E_invz[sqab > gigpara] <- as.vector(sqrt(b_z)*BesselK(sqab[sqab > gigpara], 3/2, expon.scaled = FALSE)/
                                            (sqrt(a_z[sqab > gigpara])*BesselK(sqab[sqab > gigpara], 1/2, expon.scaled = FALSE))-1/a_z[sqab > gigpara])
      E_z[sqab <= gigpara] <- as.vector(sqrt(a_z[sqab <= gigpara])*gamma(3/2)*(sqab[sqab <= gigpara]/2)^(-3/2)/
                                          ((sqrt(b_z)*gamma(1/2)*(sqab[sqab <= gigpara]/2)^(-1/2))))
      E_invz[sqab <= gigpara] <- as.vector(sqrt(b_z)*gamma(3/2)*(sqab[sqab <= gigpara]/2)^(-3/2)/
                                             (sqrt(a_z[sqab <= gigpara])*gamma(1/2)*(sqab[sqab <= gigpara]/2)^(-1/2))-1/a_z[sqab <= gigpara])
      
      # w_i
      bb <- eta2*invtau2vec*E_invsig2/2
      E_invw[1:(TForder+1)] <- (1/2+a_w)/(bb[1:(TForder+1)]+b_w)
      E_invw[(TForder+2):n] <- 1/(bb[(TForder+2):n] + E_invnu)
      E_invw[E_invw<lower] <- lower
      E_invw[E_invw>upper] <- upper
      
      # nu
      E_invnu <- 1/(2*(E_invw[(TForder+2):n] + 1/d2))
      
      # tau^2
      bbb <- sum(E_invsig2*E_invw[(TForder+2):n]*eta2[(TForder+2):n])/2 + E_invxi
      E_invtau2 <- min(max((n-TForder)/(2*bbb), 10^(-6)), 10^6)
      invtau2vec <- c(rep(1, (TForder+1)), rep(E_invtau2, (n-TForder-1)))
      
      E_invwinvtau2vec <- E_invw*invtau2vec
      E_invwinvtau2vec[E_invwinvtau2vec< 10^(-11)] <- 10^(-11)
      
      # xi
      E_invxi <- 1/(2*(E_invtau2+1))
      
    }
  }
  
  return(list(sm=th.pos, var=var.th, Esig2=E_sig2, Einvsig2=E_invsig2))
}




BQTF.VB.Lap.multi <- function(X, Y, TForder=1, pp=0.5, epsilon=10^(-3), sig.init=c(1,1)){
  
  Data <- data.frame(cbind(X,Y))
  colnames(Data) <- c("X","Y")
  Data <- Data[order(Data$X),]
  X <- Data$X
  Y <- Data$Y
  location <- unique(Data$X)
  n <- length(location)
  N <- length(Y)
  Ni <- rep(NA, n)
  y <- rep(NA, n)
  ymax <- rep(NA, n)
  num <- rep(0, n+1)
  for (i in 1:n) {
    data <- subset(Data, Data$X == location[i])
    Ni[i] <- length(data$Y)
    y[i] <- sum(data$Y)
    num[i+1] <- num[i]+Ni[i]
    ymax[i] <- max(data$Y)
  }
  
  
  ## settings
  psi <- (1-2*pp)/(pp*(1-pp))
  t2 <- 2/(pp*(1-pp))
  omega <- 1
  
  # difference operator
  D1 <- matrix(0, n-1,n)
  for (i in 1:n-1) {
    for (j in 1:n) {
      if(i == j){
        D1[i,j] <- 1
      }
      if(i+1 == j){
        D1[i,j] <- -1
      }
    }
  }
  if(TForder == 0){
    DD <- D1
    D <- rbind(c(1,rep(0, n-1)),DD)
  }
  if(TForder == 1){
    DD <- D1[1:(n-2),1:(n-1)]%*%diag(TForder/diff(location))%*%D1
    D <- rbind(t(matrix(c(1,rep(0, n-1), 
                          0,1,rep(0, n-2)), n, 2)),DD)
  }
  if(TForder == 2){
    D2 <- D1[1:(n-2),1:(n-1)]%*%diag(TForder/diff(location))%*%D1
    DD <- D1[1:(n-3),1:(n-2)]%*%diag(TForder/diff(location,lag=2))%*%D2
    D <- rbind(t(matrix(c(1,rep(0, n-1),
                          0,1,rep(0, n-2),
                          0,0, 1, rep(0, n-3)), n, 3)),DD)
  }
  
  TD <- t(D)
  invD <- solve(D)
  
  # setting
  oldtheta <- 1
  oldvar <- diag(rep(1,n))
  rr <- 1
  rrvar <- 1
  oldXX <- XX <- rep(NA, n)
  
  ## box
  th.pos <- rep(NA, n)
  var.th <- rep(NA, n)
  
  
  # initial value
  IGa <- sig.init[1]   # paramater of sig2 prior
  IGb <- sig.init[2]   # paramater of sig2 prior
  d2 <- 1  # parameter of nu prior
  lower <- 10^(-6)   # lower bound
  upper <- 10^6   # upper bound 
  zlower <- 10^(-7)   # lower bound
  zupper <- 10^7   # upper bound
  a_w <- rep(1, TForder+1)
  b_w <- rep(1, TForder+1)
  
  
  E_theta <- rep(2,n)
  E_thetai2 <- rep(4,n)
  E_invsig2 <- 1
  E_sig2 <- 1
  a_sig2 <- 1
  E_z <- rep(1,N)
  E_invz <- rep(1,N)
  a_z <- rep(1,N)
  b_z <- 1
  E_w <- rep(1,n)
  E_invw <- rep(1,n)  
  E_invnu <- 1
  E_invgamma <- 1
  E_gamma <- 1
  eta2 <- rep(1,n)
  invB <- diag(rep(1,n))
  gigpara <- 0.1
  sum.invz <- rep(1, n)
  sum.yinvz <- as.vector(sapply(split(Y*E_invz, X), sum))
  
  
  
  gigpara <- 0.08
  
  for(j in 1:100000){
    
    oldXX <- diag(oldvar)
    XX <- diag(invB/E_invsig2)
    
    rr <- max(abs(oldtheta-E_theta))
    rrvar <- max(abs(oldXX-XX))
    
    if(rr < epsilon & rrvar < epsilon){
      th.pos <- E_theta
      var.th <- diag(invB/E_invsig2)
      break
    }else{
      
      oldtheta <- E_theta
      oldvar <- invB/E_invsig2
      
      # theta
      sum.invz <- as.vector(sapply(split(E_invz, X), sum))/t2
      #sum.invz[sum.invz > upper] <- upper
      invB <- solve( diag(sum.invz) + TD%*%diag(E_invw)%*%D )
      invB <- (t(invB) + invB)/2
      C <- as.vector(sapply(split(Y*E_invz-psi, X), sum))/t2
      AA <- invB/E_invsig2 + invB%*%C%*%t(C)%*%invB
      E_theta <- invB%*%C
      E_thetai2 <- diag(AA)
      eta2 <- as.vector(diag(D%*%AA%*%TD))
      
      
      # sigma^2
      aa <- bb <- rep(0, N)
      for (i in 2:(n+1)) {
        aa[(num[i-1]+1):num[i]] <- Y[(num[i-1]+1):num[i]]^2-2*Y[(num[i-1]+1):num[i]]*E_theta[i-1]+E_thetai2[i-1]
        bb[(num[i-1]+1):num[i]] <- Y[(num[i-1]+1):num[i]]-E_theta[i-1]
      }
      a_sig2 <- sum(aa*E_invz - 2*psi*bb + psi^2*E_z)/(2*t2) + sum(eta2*E_invw)/2 + sum(E_z) + IGb
      E_invsig2 <- ((3*N+n)/2+IGa)/a_sig2
      E_sig2 <- a_sig2/((3*N+n)/2+IGa-1)
      
      
      # z
      a_z <- aa*E_invsig2/t2
      b_z <- (psi^2*E_invsig2/t2+2*E_invsig2)
      sqab <- sqrt(a_z*b_z)
      E_z[sqab > gigpara] <- as.vector(sqrt(a_z[sqab > gigpara])*BesselK(sqab[sqab > gigpara], 3/2, expon.scaled = FALSE)/
                                         (sqrt(b_z)*BesselK(sqab[sqab > gigpara], 1/2, expon.scaled = FALSE)))
      E_invz[sqab > gigpara] <- as.vector(sqrt(b_z)*BesselK(sqab[sqab > gigpara], 3/2, expon.scaled = FALSE)/
                                            (sqrt(a_z[sqab > gigpara])*BesselK(sqab[sqab > gigpara], 1/2, expon.scaled = FALSE))-1/a_z[sqab > gigpara])
      E_z[sqab <= gigpara] <- as.vector(sqrt(a_z[sqab <= gigpara])*gamma(3/2)*(sqab[sqab <= gigpara]/2)^(-3/2)/
                                          ((sqrt(b_z)*gamma(1/2)*(sqab[sqab <= gigpara]/2)^(-1/2))))
      E_invz[sqab <= gigpara] <- as.vector(sqrt(b_z)*gamma(3/2)*(sqab[sqab <= gigpara]/2)^(-3/2)/
                                             (sqrt(a_z[sqab <= gigpara])*gamma(1/2)*(sqab[sqab <= gigpara]/2)^(-1/2))-1/a_z[sqab <= gigpara])
      
      
      # w_i
      A_w <- E_gamma
      B_w <- E_invsig2*eta2
      sqAB <- sqrt(A_w*B_w)
      E_invw[1:(TForder+1)] <- (1/2+a_w)/(B_w[1:(TForder+1)]/2+b_w)
      E_w[1:(TForder+1)] <- (B_w[1:(TForder+1)]/2+b_w)/(a_w-1/2)
      E_invw[(TForder+2):n] <- as.vector((BesselK(sqAB[(TForder+2):n], 3/2, expon.scaled = FALSE)/BesselK(sqAB[(TForder+2):n], 1/2, expon.scaled = FALSE))*
                                           (A_w/B_w[(TForder+2):n])^(1/2) - 1/B_w[(TForder+2):n] )
      E_w[(TForder+2):n] <- as.vector((BesselK(sqAB[(TForder+2):n], 3/2, expon.scaled = FALSE)/BesselK(sqAB[(TForder+2):n], 1/2, expon.scaled = FALSE))*
                                        (B_w[(TForder+2):n]/A_w)^(1/2))
      E_invw[E_invw<lower] <- lower
      E_invw[E_invw>upper] <- upper
      
      
      # gam2
      a_gamma <- 2*E_invnu
      b_gamma <- sum(E_w[(TForder+2):n])
      sqab <- sqrt(a_gamma*b_gamma)
      if(sqab > gigpara){
        E_gamma <- sqrt(a_gamma)*BesselK(sqab, n-TForder-1/2, expon.scaled = FALSE)/(sqrt(b_gamma)*BesselK(sqab, n-TForder-3/2, expon.scaled = FALSE))
        E_invgamma <- sqrt(b_gamma)*BesselK(sqab, n-TForder-1/2, expon.scaled = FALSE)/(sqrt(a_gamma)*BesselK(sqab, n-TForder-3/2, expon.scaled = FALSE)) - 2*(n-TForder-3/2)/a_gamma
      }else{
        E_gamma <- sqrt(a_gamma)*gamma(n-TForder-1/2)*(sqab/2)^(-n+TForder+1/2)/(sqrt(b_gamma)*gamma(n-TForder-3/2)*(sqab/2)^(-n+TForder+3/2))
        E_invgamma <- sqrt(b_gamma)*gamma(n-TForder-1/2)*(sqab/2)^(-n+TForder+1/2)/(sqrt(a_gamma)*gamma(n-TForder-3/2)*(sqab/2)^(-n+TForder+3/2)) - 2*(n-TForder-3/2)/a_gamma
        
      }
      
      
      # nu
      E_invnu <- 1/(2*(E_invgamma + 1/d2))
      
    }
    #if(round(j/200)==(j/200)){ print(j) }
  }
  
  return(list(sm=th.pos, var=var.th, E_sig2=E_sig2, E_invsig2=E_invsig2,
              E_invw=E_invw, E_w=E_w, E_invnu=E_invnu, E_z=E_z, E_invz=E_invz,
              E_invnu=E_invnu, E_invgamma=E_invgamma, E_gamma=E_gamma))
}






BQTF.CVB.multi <- function(X, Y, TForder=1, pp=0.5, epsilon=c(10^(-3),10^(-2)), prior="HS", B.size=100, sig.init=c(0.1,0.1), type=1){
  
  eps1 <- epsilon[1]
  eps2 <- epsilon[2]
  Data <- data.frame(cbind(X,Y))
  colnames(Data) <- c("X","Y")
  Data <- Data[order(Data$X),]
  X <- Data$X
  Y <- Data$Y
  location <- unique(Data$X)
  n <- length(location)
  N <- length(Y)
  Ni <- rep(NA, n)
  y <- rep(NA, n)
  ymax <- rep(NA, n)
  num <- rep(0, n+1)
  for (i in 1:n) {
    data <- subset(Data, Data$X == location[i])
    Ni[i] <- length(data$Y)
    y[i] <- sum(data$Y)
    num[i+1] <- num[i]+Ni[i]
    ymax[i] <- max(data$Y)
  }
  
  if(prior == "HS"){
    fit <- BQTF.VB.HS.multi(X, Y, TForder=TForder, pp=pp, epsilon= eps1, sig.init=sig.init)
    var <- fit$var
    sm <- fit$sm
    if(pp !=0.5){
      fit05 <- BQTF.VB.HS.multi(X, Y, TForder=TForder, pp=0.5, epsilon= eps1, sig.init=sig.init)
      sm05 <- fit05$sm
    }else{
      sm05 <- sm
    }
  }
  if(prior == "Lap"){
    fit <- BQTF.VB.Lap.multi(X, Y,  TForder=TForder, pp=pp, epsilon=eps1, sig.init=sig.init)
    var <- fit$var
    sm <- fit$sm
    if(pp !=0.5){
      fit05 <- BQTF.VB.Lap.multi(X, Y,  TForder=TForder, pp=0.5, epsilon=eps1, sig.init=sig.init)
      sm05 <- fit05$sm
    }else{
      sm05 <- sm
    }
  }
  
  ## bootstrap samples (sample size: B.size)
  zansa <- rep(NA, N)
  for (i in 1:n) {
    data <- subset(Data, Data$X == location[i])
    zansa[(num[i]+1):num[i+1]] <- data$Y - sm05[i]
  }
  Bsample <- Bdata <- matrix(NA, N,  B.size)
  for (bb in 1:B.size) {
    Bsample[,bb] <- sample(zansa, N, replace = TRUE)
  }
  for (i in 1:n) {
    Bdata[(num[i]+1):num[i+1],] <- sm05[i] + Bsample[(num[i]+1):num[i+1],]
  }
  
  sample <- Result1var <- matrix(NA, B.size, n)
  
  ## bootstrap sample
  for (bb in 1:B.size) {
    
    if(prior == "HS"){
      Result1 <- BQTF.VB.HS.multi(X, Bdata[,bb], TForder=TForder, pp=pp, epsilon= eps2, sig.init=sig.init)
    }
    if(prior == "Lap"){
      Result1 <- BQTF.VB.Lap.multi(X, Bdata[,bb], TForder=TForder, pp=pp, epsilon=eps2, sig.init=sig.init)
    }
    
    sample[bb,] <- as.vector(Result1$sm)
    Result1var[bb,] <- as.vector(Result1$var)
    
    if(round(bb/20)==(bb/20)){ print(bb) }
  }
  
  
  ## tuning of lambda
  if(type==1){
    CP <- lambda <- c()
    selectedvar  <- 100000001*var 
    selectedlambda <- 100000001
    for(j in 1:100000000){
      lambda[j] <- 1+(0.1)*(j-1)
      CI <- matrix(NA, n,2)
      CI[,1] <- qnorm(0.025, sm, sqrt(lambda[j]*var))
      CI[,2] <- qnorm(0.975, sm, sqrt(lambda[j]*var))
      count <- rep(0, n)
      for(i in 1:n){
        count[i] <- sum( CI[i,1] <= sample[,i] & sample[,i] <= CI[i,2]) 
      }
      CP[j] <- sum(count)/(B.size*n)
      
      if(CP[j] >= 0.95){
        if(j == 1){
          selectedlambda <- lambda[j]
        }else{
          selectedlambda <- (c(lambda[j],lambda[j-1])[order(abs(CP[j]-0.95), abs(CP[j-1]-0.95))])[1]
        }
        selectedvar <- selectedlambda*var
        break
      }
    }
  }
  if(type==2){
    CP <- lambda <- c()
    selectedvar <- rep(1, n)
    selectedlambda <- rep(1, n)
    for (i in 1:n) {
      for(j in 1:100000000){
        lambda[j] <- 1+(0.1)*(j-1)
        CI <- matrix(NA, n,2)
        CI1 <- qnorm(0.025, sm[i], sqrt(lambda[j]*var[i]))
        CI2 <- qnorm(0.975, sm[i], sqrt(lambda[j]*var[i]))
        count <- sum( CI1 <= sample[,i] & sample[,i] <= CI2 ) 
        CP[i] <- count/B.size
        if(CP[i] >= 0.95){
          if(j == 1){
            selectedlambda[i] <- lambda[j]
          }else{
            selectedlambda[i] <- (c(lambda[j],lambda[j-1])[order(abs(CP[j]-0.95), abs(CP[j-1]-0.95))])[1]
          }
          break
        }
      }
      selectedvar[i] <- selectedlambda[i]*var[i]
    }
  }
  
  
  return(list(sm=sm, var=selectedvar, selectedlambda=selectedlambda,
              origialvar=var, th.pos.sample=sample, th.pos.var=Result1var))
}

