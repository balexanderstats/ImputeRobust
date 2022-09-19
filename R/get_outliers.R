#these functions are taken from the extremevalues package on CRAN
#extreme values package author: Mark van der Loo <mark.vanderloo@gmail.com>


fitLognormal <- function(y, p)
{
  if ( !is.vector(y) )
    stop("First argument is not of type vector")
  if ( sum(y<=0) > 0 )
    stop("First argument contains nonpositive values")
  if ( !is.vector(p))
    stop("First argument is not of type vector")
  if ( sum(p<=0) > 0 | sum(p>=1) >0 )
    stop("Second argument contains values out of range (0,1)")
  if (length(y) != length(p))
    stop("First and second argument have different length");

  N <- length(y);
  lnY <- as.matrix(log(y),nrow=N)
  p <- as.matrix(p,nrow=N)


  A <- matrix(0,nrow=N,ncol=2)
  A[,1] <- 1+double(N);
  A[,2] <- sqrt(2)*invErf(2*p-1)
  param <- solve(t(A) %*% A)%*%t(A)%*%lnY
  r2 <- 1 - var(exp(A%*%param) - y)/var(y);

  return(list(mu=param[1], sigma=param[2], R2=r2));
}
fitExponential <- function(y,p)
{
  if ( !is.vector(y) )
    stop("First argument is not of type vector")
  if ( sum(y<=0) > 0 )
    stop("First argument contains nonpositive values")
  if ( !is.vector(p))
    stop("First argument is not of type vector")
  if ( sum(p<=0) > 0 | sum(p>=1) >0 )
    stop("Second argument contains values out of range (0,1)")
  if (length(y) != length(p))
    stop("First and second argument have different length");


  #   Lambda <- -sum(log(1-p))/sum(y);
  Lambda <- -sum(log(1-p)^2)/sum(y * log(1-p))
  r2 <- 1 - var(y - (-log(1-p)/Lambda) )/var(y);

  return(list(lambda=Lambda, R2=r2));


}

fitWeibull <- function(y, p)
{
  if ( !is.vector(y) )
    stop("First argument is not of type vector")
  if ( sum(y<=0) > 0 )
    stop("First argument contains nonpositive values")
  if ( !is.vector(p))
    stop("First argument is not of type vector")
  if ( sum(p<=0) > 0 | sum(p>=1) >0 )
    stop("Second argument contains values out of range (0,1)")
  if (length(y) != length(p))
    stop("First and second argument have different length");

  N <- length(y);
  lnY <- as.matrix(log(y),nrow=N)
  p <- as.matrix(p,nrow=N)


  A <- matrix(0,nrow=N,ncol=2)
  A[,1] <- 1+double(N);
  A[,2] <- log(log(1/(1-p)))
  param <- solve(t(A) %*% A)%*%t(A)%*%lnY
  r2 <- 1 - var(exp(A%*%param) - y)/var(y);

  return(list(k=1/param[2], lambda=exp(param[1]), R2=r2));
}


fitPareto <- function(y,p)
{
  if ( !is.vector(y) )
    stop("First argument is not of type vector")
  if ( sum(y<=0) > 0 )
    stop("First argument contains nonpositive values")
  if ( !is.vector(p))
    stop("First argument is not of type vector")
  if ( sum(p<=0) > 0 | sum(p>=1) >0 )
    stop("Second argument contains values out of range (0,1)")
  if (length(y) != length(p))
    stop("First and second argument have different length");

  N <- length(y);
  lnY <- as.matrix(log(y),nrow=N)
  p <- as.matrix(p,nrow=N)

  A <- matrix(0,nrow=N,ncol=2)
  A[,1] <- 1 + double(N);
  A[,2] <- log(1-p);
  param <- solve(t(A) %*% A)%*%t(A)%*%lnY
  r2 <- 1 - var(exp(A%*%param) - y)/var(y);

  return(list(ym=exp(param[1]), alpha=-1/param[2], R2=r2));
}

fitNormal <- function(y, p)
{
  if ( !is.vector(y) )
    stop("First argument is not of type vector")
  if ( !is.vector(p))
    stop("First argument is not of type vector")
  if ( sum(p<=0) > 0 | sum(p>=1) >0 )
    stop("Second argument contains values out of range (0,1)")
  if (length(y) != length(p))
    stop("First and second argument have different length");

  N <- length(y);
  Y <- as.matrix(y,nrow=N)
  p <- as.matrix(p,nrow=N)


  A <- matrix(0,nrow=N,ncol=2)
  A[,1] <- 1+double(N);
  A[,2] <- sqrt(2)*invErf(2*p-1)
  param <- solve(t(A) %*% A) %*% t(A) %*% Y
  r2 <- 1 - var(A%*%param - y)/var(y);
  return(list(mu=param[1], sigma=param[2], R2=r2));
}

getNormalLimit <- function(y, p, N, rho)
{
  param <- fitNormal(y,p)
  ell <- c(Left=-Inf, Right=Inf)
  if ( !is.na(rho[1]) )
    ell[1] <- sqrt(2)*param$sigma*invErf(2*rho[1]/N-1)+param$mu
  if ( !is.na(rho[2]) )
    ell[2] <- sqrt(2)*param$sigma*invErf(1-2*rho[2]/N)+param$mu
  return(list(mu=param$mu,
              sigma=param$sigma,
              nFit=length(y),
              R2=param$R2,
              limit=ell)
  )
}


getExponentialLimit <- function(y, p, N, rho)
{
  param <- fitExponential(y, p)
  ell <- c(Left=0, Right=Inf)
  if ( !is.na(rho[1]) )
    ell[1] <- -log(1-rho[1]/N)/param$lambda
  if ( !is.na(rho[2]) )
    ell[2] <- log(N/rho[2])/param$lambda

  return(list(lambda=param$lambda,
              R2=param$R2,
              nFit=length(y),
              limit=ell)
  )
}





# wrapper function for outlier detection methods.
# 23.12.2009 version 1, mvdl
getOutliers <- function(y, method="I",  ...)
{
  if ( !(method %in% c("I","II") ) )
    stop("method not recognized (choose I or II)")
  out <- switch( method,
                 I = getOutliersI(y, ...)#,
                 #II = getOutliersII(y, ...)
  )
  return(out)
}
invErf <- function(x)
{
  if ( sum(x >= 1) > 0  | sum(x <= -1) > 0 )
    stop("Argument must be between -1 and 1")

  return(qnorm((1+x)/2)/sqrt(2));

}

# 19.10.2009, changed p-estimator
# 24.11.2009, added Weibull distribution
# 22.12.2009, added left limit, Changed rho default, added input checks.
# 08.01.2010, switched to Makkonen's equation for plot positions.
getOutliersI <- function(y, rho=c(1,1), FLim=c(0.1,0.9), distribution="normal")
{

  if ( !is.vector(y) )
    stop("First argument is not of type vector")
  if ( sum(y < 0) > 0 & !(distribution == "normal") )
    stop("First argument contains nonpositive values")
  if ( sum( rho <= 0, na.rm=TRUE ) > 0 )
    stop("Values of rho must be positive")
  if ( FLim[2] <= FLim[1] | sum( FLim < 0 | FLim > 1) >0 )
    stop("Invalid range in FLim: 0<=FLim[1]<FLim[2]<=1")
  if ( ! distribution %in% c("lognormal", "pareto", "exponential", "weibull", "normal") )
    stop("Invalid distribution (lognormal, pareto, exponential, weibull, normal).")

  Y <- y;

  y <- sort(y);
  N <- length(y)
  P <- (1:N)/(N+1)
  Lambda <- P >= FLim[1] & P<=FLim[2]

  y <- y[Lambda];
  p <- P[Lambda];
  out <- switch(distribution,
                lognormal = getLognormalLimit(y, p, N, rho),
                pareto = getParetoLimit(y, p, N, rho),
                exponential = getExponentialLimit(y, p, N, rho),
                weibull = getWeibullLimit(y, p, N, rho),
                normal = getNormalLimit(y, p, N, rho)
  )

  out$method <- "Method I"
  out$distribution=distribution
  out$iRight = which( Y > out$limit[2] )
  out$iLeft = which( Y < out$limit[1] )
  out$nOut = c(Left=length(out$iLeft), Right=length(out$iRight))
  out$yMin <- y[1]
  out$yMax <- tail(y,1)
  out$rho = c(Left=rho[1], Right=rho[2])
  out$Fmin = FLim[1]
  out$Fmax = FLim[2]

  return(out);
}


getParetoLimit <- function(y, p, N, rho)
{
  param <- fitPareto(y,p)
  ell <- c(Left=-Inf, Right=Inf)
  if ( !is.na(rho[1]) )
    ell[1] <- param$ym*(1-rho[2]/N)^{-1/param$alpha}
  if ( !is.na(rho[2]) )
    ell[2] <- param$ym*(N/rho[2])^{1/param$alpha}

  return(list(ym=param$ym,
              alpha=param$alpha,
              nFit=length(y),
              R2=param$R2,
              limit=ell)
  )
}




getWeibullLimit <- function(y, p, N, rho)
{
  param <- fitWeibull(y,p)
  ell <- c(Left=0, Right=Inf)
  if ( !is.na(rho[1]) )
    ell[1] <- param$lambda*(-log(1-rho[1]/N))^(1/param$k)
  if ( !is.na(rho[2]) )
    ell[2] <- param$lambda*log(N/rho[2])^(1/param$k)
  return(list(k=param$k,
              lambda=param$lambda,
              nFit=length(y),
              R2=param$R2,
              limit=ell)
  )
}




getLognormalLimit <- function(y, p, N, rho)
{
  param <- fitLognormal(y,p)
  ell <- c(Left=0, Right=Inf)
  if ( !is.na(rho[1]) )
    ell[1] <- exp(sqrt(2)*param$sigma*invErf(-(1-2*rho[1]/N))+param$mu)
  if ( !is.na(rho[2]) )
    ell[2] <- exp(sqrt(2)*param$sigma*invErf(1-2*rho[2]/N)+param$mu)
  return(list(mu=param$mu,
              sigma=param$sigma,
              nFit=length(y),
              R2=param$R2,
              limit=ell)
  )
}


