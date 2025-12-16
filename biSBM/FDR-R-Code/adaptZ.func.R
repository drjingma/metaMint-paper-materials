#' @param model implementation models. Could be one of the following
#' \describe{
#'   \item{\code{Gauss}}{all Gaussian parameters of the null distribution are unknown and estimated from data.}
#'   \item{\code{Gauss01}}{Compared to \code{Gauss}, the null distribution is set to N(0,1).}
#' }
# Heavy step: run once -----------------------------------------------
adaptZ_prepare <- function(zv, gamma, model = c("Gauss", "Gauss01"), density_n = 1000L) {
  model <- match.arg(model)
  
  # 1) kernel density on a wide support, then interpolate at zv
  zv.ds.obj <- density(zv, from = min(zv) - 10, to = max(zv) + 10, n = density_n)
  zv.ds     <- lin.itp(zv, zv.ds.obj$x, zv.ds.obj$y)   # interpolate density at each zv
  
  # 2) null params
  if (model == "Gauss") {
    zs <- EstNull.func(zv, gamma)  # returns list(mu=..., s=...)
    mu <- zs$mu; s <- zs$s
  } else { # "Gauss01"
    mu <- 0; s <- 1
  }
  
  # 3) mixture weight and local fdr
  zv.p0   <- 1 - epsest.func(zv, mu, s)               # π0(z) = 1 - ε(z)
  zv.lfdr <- zv.p0 * dnorm(zv, mean = mu, sd = s) / zv.ds
  
  # clip numerical extremes just in case (optional)
  zv.lfdr <- pmin(pmax(zv.lfdr, 0), 1)
  
  # return everything needed for thresholding
  list(
    zv      = zv,
    mu      = mu,
    s       = s,
    dens    = zv.ds,
    p0      = zv.p0,
    lfdr    = zv.lfdr
  )
}

# Light step: apply any number of q's --------------------------------
adaptZ_apply <- function(prep, q) {
  # adpt.cutz() is your existing thresholding routine on lfdr
  y <- adpt.cutz(prep$lfdr, q)
  # keep original return structure
  y
}

# Convenience: vector of q's -> list of results ----------------------
adaptZ_apply_many <- function(prep, q_vec) {
  setNames(lapply(q_vec, function(q) adaptZ_apply(prep, q)), paste0("q=", q_vec))
}

adaptZ.func<-function(zv, q, gamma, model){
  # the input
    # zv is the z-values transformed from m tests
    # q is the desired FDR level
  # the output is a list with
    # the first element (st.lfdr) the sorted local fdr values
    # the second element (k) the number of hypotheses to be rejected
    # the third element (lfdrk) the threshold for the local fdr values
    # the fourth element (reject) the set of indices of the rejected hypotheses
    # the fifth element (accept) the set of indices of the accepted hypotheses
   ## the estimates for the local fdr statistics

  # density estimates
  zv.ds<-density(zv, from=min(zv)-10, to=max(zv)+10, n=1000)
  # linear interpolation
  zv.ds<-lin.itp(zv, zv.ds$x, zv.ds$y)
  if (model=='Gauss'){
    # estimating the null distribution
    # estimates for the null density are sensitive to the gamma parameter
    zv.MuSigma<-EstNull.func(zv,gamma)
    mu<-zv.MuSigma$mu
    s<-zv.MuSigma$s
  } else if (model=='Gauss01'){
    mu<-0; s<-1
  }
  zv.p0<-1-epsest.func(zv, mu, s)
  zv.lfdr<-zv.p0*dnorm(zv, mu, s)/zv.ds
  y<-adpt.cutz(zv.lfdr, q)
  return (y)
}

epsest.func <- function(x,u,sigma)
{
  # x is a vector
  # u is the mean
  # sigma is the standard deviation

  z  = (x - u)/sigma
  xi = c(0:100)/100
  tmax=sqrt(log(length(x)))
  tt=seq(0,tmax,0.1)

  epsest=NULL

  for (j in 1:length(tt)) {

    t  = tt[j]
    f  = t*xi
    f  = exp(f^2/2)
    w  = (1 - abs(xi))
    co  = 0*xi

    for (i in 1:101) {
      co[i] = mean(cos(t*xi[i]*z));
    }
    epshat = 1 - sum(w*f*co)/sum(w)
    epsest = c(epsest,epshat)
  }
  return(epsest=max(epsest))
}

EstNull.func<-function (x,gamma=0.1)
{
  # x is a vector of z-values
  # gamma is a parameter, default is 0.1
  # output the estimated mean and standard deviation

  n = length(x)
  t = c(1:1000)/200

  gan    = n^(-gamma)
  that   = 0
  shat   = 0
  uhat   = 0
  epshat = 0

  phiplus   = rep(1,1000)
  phiminus  = rep(1,1000)
  dphiplus  = rep(1,1000)
  dphiminus = rep(1,1000)
  phi       = rep(1,1000)
  dphi      = rep(1,1000)

  for (i in 1:1000) {
    s = t[i]
    phiplus[i]   = mean(cos(s*x))
    phiminus[i]  = mean(sin(s*x))
    dphiplus[i]  = -mean(x*sin(s*x))
    dphiminus[i] = mean(x*cos(s*x))
    phi[i]       = sqrt(phiplus[i]^2 + phiminus[i]^2)
  }

  ind = min(c(1:1000)[(phi - gan) <= 0])
  tt = t[ind]
  a  = phiplus[ind]
  b  = phiminus[ind]
  da = dphiplus[ind]
  db = dphiminus[ind]
  c  = phi[ind]

  that   = tt
  shat   = -(a*da + b*db)/(tt*c*c)
  shat   = sqrt(shat)
  uhat   = -(da*b - db*a)/(c*c)
  epshat = 1 - c*exp((tt*shat)^2/2)

  return(musigma=list(mu=uhat,s=shat))
}

lin.itp<-function(x, X, Y){
  ## x: the coordinates of points where the density needs to be interpolated
  ## X: the coordinates of the estimated densities
  ## Y: the values of the estimated densities
  ## the output is the interpolated densities
  x.N<-length(x)
  X.N<-length(X)
  y<-rep(0, x.N)
  for (k in 1:x.N){
    i<-max(which((x[k]-X)>=0))
    if (i<X.N)
      y[k]<-Y[i]+(Y[i+1]-Y[i])/(X[i+1]-X[i])*(x[k]-X[i])
    else
      y[k]<-Y[i]
  }
  return(y)
}

adpt.cutz<-function(lfdr, q)
{
  # the input
  # lfdr the vector of local fdr statistics
  # q the desired FDR level
  # the output is a list with
  # the first element (st.lfdr) the sorted local fdr values
  # the second element (k) the number of hypotheses to be rejected
  # the third element (lfdrk) the threshold for the local fdr values
  # the fourth element (reject) the set of indices of the rejected hypotheses
  # the fifth element (accept) the set of indices of the accepted hypotheses

  m=length(lfdr)
  st.lfdr<-sort(lfdr)
  k=1
  while(k<m && (1/k)*sum(st.lfdr[1:k])<q){
    k=k+1
  }
  k<-k-1
  lfdrk<-st.lfdr[k]
  reject<-which(lfdr<=lfdrk)
  accept<-which(lfdr>lfdrk)
  y<-list(sf=st.lfdr, nr=k, thr=lfdrk, re=reject, ac=accept)
  return (y)
}
