# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'


calI.cont <- function(dat1){
  weight <- 1/dat1$s2
  theta_f <- weighted.mean(dat1$y,weight)
  Q <- sum(weight*(dat1$y-theta_f)^2)
  H2 <- Q/(nrow(dat1)-1)
  I2 <- (Q>(nrow(dat1)-1))*(1-1/H2)
  return(I2)
}


meta.reml <- function(dat){
  # obtain reml estimator of tau
  return(
    tryCatch(
      {
        res3 <- metafor::rma(yi=y,vi=s2,data=dat,method='REML',control = list(stepadj=0.5,maxiter=1000))
        return(res3)
      },
      error=function(e) {
        dat$study.name <- c(1:nrow(dat))
        res3 <- metafor::rma.mv(y, s2, random = ~ 1 | study.name, data=dat, control=list(optimizer="optim", optmethod="Nelder-Mead"))
        res3$tau2 <- res3$sigma2
        return(res3)}
    ))
}

#' @export
calR <- function(dat1,m){
  # obtain Rm statistic
  lev_set <- combn(nrow(dat1),m)
  ncombn <- ncol(lev_set)
  Rval_ls <- sapply(c(1:ncombn),function(c){
    dat2 <- dat1[-lev_set[,c],]
    weight <- 1/dat2$s2
    theta_f <- weighted.mean(dat2$y,weight)
    res <- meta.reml(dat2)
    tau2 <- res$tau2
    var.tr <- res$se^2
    theta_r <- res$beta
    dat3 <- dat1[lev_set[,c],]
    d_i <- dat3$y-theta_r
    var_di <- dat3$s2+tau2+var.tr
    ri <- d_i/sqrt(var_di)
    R_value <- sum(ri^2)/m
    return(R_value)
  })
  return(list(maxR=max(Rval_ls),Rvals=Rval_ls))
}

bootdat <- function(dat1,SEED,tau2){
  set.seed(SEED)
  dat1$s2 <- sample(dat1$s2,nrow(dat1),replace=TRUE)
  w_as <- 1/(dat1$s2+tau2)
  theta_r <- weighted.mean(dat1$y,w_as)
  dat1.boot <- data.frame(y=rnorm(nrow(dat1), theta_r, sqrt(dat1$s2 + tau2)), s2=dat1$s2)
  return(dat1.boot)
}

#' @export
calempF <- function(dat1,m){
  set.seed(9876)
  tau2 <- meta.reml(dat1)$tau2
  ## 1000 bootstrap runs
  seed.ls <- sample(1:200000,1000)
  Rmls <- parallel::mclapply(seed.ls,function(SEED){
    dat1.boot <- bootdat(dat1,SEED,tau2)
    ## leave one study out
    maxR <- calR(dat1.boot,m)
    Rm <- maxR$maxR
    return(Rm)
  },mc.cores=6)
  ## use bootstrap to estimate the sample variance
  Rval2 <- data.frame(do.call(cbind,Rmls))
  Rval2 <- as.vector(Rval2)
  return(Rval2)
}

#' @export
Rm.func.iterative.boot <- function(dat1,m=1){
  # function to detect and identify non-replicable studies for bootstrap method
  Rls <- calR(dat1,m)
  Rls$empdist <- calempF(dat1,m)
  Pval <- (sum(Rls$empdist>Rls$maxR)+1)/(length(Rls$empdist)+1)
  Pval_leave1 <- Pval
  dat1.copy <- dat1
  out.n <- 1
  while(Pval<0.05 & out.n < (nrow(dat1.copy))){
    top2.ind <- order(Rls$Rvals,decreasing = T)[1:2]
    R.top2 <- sapply(top2.ind,function(x){
      dat1.sub <- dat1[-x,]
      Rls <- calR(dat1,1)
      return(Rls$maxR)
    })
    sm.ind <- which.min(R.top2)
    dat1.sub <- dat1[-top2.ind[sm.ind],]
    Rls <- calR(dat1.sub,1)
    Pval <- (sum(Rls$empdist>Rls$maxR)+1)/(length(Rls$empdist)+1)
    out.n <- out.n+1*(Pval<0.05)
    dat1 <- dat1.sub
  }
  out.stdies <- which(!dat1.copy$y %in% dat1$y)
  return(list(stat=c(Rls$maxR,Pval),P_1=Pval_leave1,out_studies=out.stdies))
}

#' @export
Rm.func.iterative <- function(dat1,m=1){
  # function to detect and identify non-replicable study under Gumbel approximation
  Rls <- calR(dat1,m)
  lev_set <- combn(nrow(dat1),m)
  ncombn <- ncol(lev_set)
  cn <- 2
  dn=2*(log(ncombn)+(m/2-1)*log(log(ncombn))-log(gamma(m/2)))
  Pval <- evd::pgumbel((1/cn)*(m*Rls$maxR-dn),loc=0,lower.tail = F)
  Pval_leave1 <- Pval
  dat1.copy <- dat1
  out.n <- 1
  # while(Pval<0.05 & out.n < min(4,floor(nrow(dat1)/2))){
  while(Pval<0.05 & out.n < (nrow(dat1.copy))){
    top2.ind <- order(Rls$Rvals,decreasing = T)[1:2]
    R.top2 <- sapply(top2.ind,function(x){
      dat1.sub <- dat1[-x,]
      Rls <- calR(dat1,1)
      return(Rls$maxR)
    })
    sm.ind <- which.min(R.top2)
    dat1.sub <- dat1[-top2.ind[sm.ind],]
    Rls <- calR(dat1.sub,1)
    lev_set <- combn(nrow(dat1.sub),m)
    ncombn <- ncol(lev_set)
    cn <- 2
    dn=2*(log(ncombn)+(m/2-1)*log(log(ncombn))-log(gamma(m/2)))
    Pval <- evd::pgumbel((1/cn)*(m*Rls$maxR-dn),loc=0,lower.tail = F)
    out.n <- out.n+1*(Pval<0.05)
    dat1 <- dat1.sub
  }
  out.stdies <- which(!dat1.copy$y %in% dat1$y)
  return(list(stat=c(Rls$maxR,Pval),P_1=Pval_leave1,out_studies=out.stdies))
}
