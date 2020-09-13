titecrm <- function(prior, target, tox, level, n=length(level),
                    weights=NULL, followup=NULL, entry=NULL, exit=NULL,
                    obswin=NULL, scheme="linear", conf.level=0.90,
                    dosename=NULL, include=1:n, pid=1:n, method="bayes",model="empiric",var.est=TRUE,
                    scale=sqrt(1.34), intcpt=3, model.detail=TRUE, patient.detail=TRUE, tite=TRUE) {
  if (is.null(weights)) {
    if (is.null(followup)) { followup <- exit-entry; }
    if (scheme=="linear") { weights <- followup/obswin; }
    else if (scheme=="adaptive") {
      support <- sort(followup[tox==1])
      z <- length(support)
      if (z) {
        for (i in 1:n) {
          m <- length(support[support<=followup[i]])
          if (!m) weights[i] <- followup[i] / support[1] / (z+1)
          else if (m==z) weights[i] <- (z + (followup[i]-support[z])/(obswin-support[z]))/(z+1)
          else weights[i] <- (m + (followup[i]-support[m])/(support[m+1]-support[m]))/(z+1)
        }
      }
      else { weights <- followup/obswin; }
    }
    else { stop(" Weighting scheme undefined!"); }
    weights <- pmin(weights, 1)
  }

  if (any(weights>1) | any(weights<0)) stop(" Weights have to be between 0 and 1!")

  if (is.null(pid)) {
    if (! (length(tox)==length(level) & length(tox)==length(weights)))
      stop(" tox, level, and weights are of different lengths!")
  }
  else {
    if (! (length(tox)==length(level) & length(tox)==length(weights) & length(tox)==length(pid)) )
      stop(" pid, tox, level, and weights are of different lengths!")
  }
  weights[tox==1] <- 1

  y1p <- tox[include]
  w1p <- weights[include]
  if (model=="empiric") {
    dosescaled <- prior

    x1p <- prior[level[include]]
    if (method=="mle") {
      if (sum(y1p)==0 | sum(y1p)==length(y1p)) stop(" mle does not exist!")
      est <- optimize(lcrm,c(-10,10),x1p,y1p,w1p,tol=0.0001,maximum=TRUE)$max
      if (var.est) { e2 <- integrate(crmht2,-10,10,x1p,y1p,w1p,500,abs.tol=0)[[1]] / integrate(crmh,-10,10,x1p,y1p,w1p,500,abs.tol=0)[[1]]; }
    }
    else if (method=="bayes") {
      den <- integrate(crmh,-Inf,Inf,x1p,y1p,w1p,scale,abs.tol=0)[[1]]
      est <- integrate(crmht,-10,10,x1p,y1p,w1p,scale,abs.tol=0)[[1]] / den
      if (var.est) { e2 <- integrate(crmht2,-10,10,x1p,y1p,w1p,scale,abs.tol=0)[[1]] / den; }
    }
    else { stop(" unknown estimation method"); }
    ptox <- prior^exp(est)
    if (var.est) {
      post.var <- e2-est^2
      crit <- qnorm(0.5+conf.level/2)
      lb <- est - crit*sqrt(post.var)
      ub <- est + crit*sqrt(post.var)
      ptoxL <- prior^exp(ub)
      ptoxU <- prior^exp(lb)
    }
  }
  else if (model=="logistic") {
    dosescaled <- log(prior/(1-prior)) - intcpt
    if (!all(dosescaled<0)) {
      stop( "intercept parameter in logit model is too small: scaled doses > 0!")
    }
    #    LB <- log(  (log((1+target)/(1-target)) - intcpt)/dosescaled[1] ) - 3
    #    UB <- log(  (log((target/2)/(1-target/2)) - intcpt)/dosescaled[length(prior)] ) + 3

    x1p <- dosescaled[level[include]]

    if (method=="mle") {
      if (sum(y1p)==0 | sum(y1p)==length(y1p)) stop(" mle does not exist!")
      est <- optimize(lcrmlgt,c(-10,10),x1p,y1p,w1p,intcpt,tol=0.0001,maximum=TRUE)$max
      if (var.est) { e2 <- integrate(crmht2lgt,-10,10,x1p,y1p,w1p,500,intcpt,abs.tol=0)[[1]] / integrate(crmhlgt,-10,10,x1p,y1p,w1p,500,abs.tol=0)[[1]]; }
    }
    else if (method=="bayes") {
      den <- integrate(crmhlgt,-Inf,Inf,x1p,y1p,w1p,scale,intcpt,abs.tol=0)[[1]]
      est <- integrate(crmhtlgt,-10,10,x1p,y1p,w1p,scale,intcpt,abs.tol=0)[[1]] / den
      if (var.est) { e2 <- integrate(crmht2lgt,-10,10,x1p,y1p,w1p,scale,intcpt,abs.tol=0)[[1]] / den; }
    }
    else { stop(" unknown estimation method"); }
    #    est <- min(UB, max(LB, est))
    ptox <- (1 + exp(-intcpt-exp(est)*dosescaled))^{-1}
    if (var.est) {
      post.var <- e2-est^2
      crit <- qnorm(0.5+conf.level/2)
      lb <- est - crit*sqrt(post.var)
      ub <- est + crit*sqrt(post.var)
      ptoxL <- (1 + exp(-intcpt-exp(ub)*dosescaled))^{-1}
      ptoxU <- (1 + exp(-intcpt-exp(lb)*dosescaled))^{-1}
    }
  }
  else { stop(" model specified not available."); }

  if (all(ptox<=target)) { rec <- length(prior); }
  else if (all(ptox>=target)) { rec <- 1; }
  else { rec <- order(abs(ptox-target))[1]; }
  if (!var.est) { post.var <- ptoxL <- ptoxU <- NA; }
  foo <- list(prior=prior, target=target, tox=tox, level=level,
              dosename=dosename, subset=pid[include], estimate=est,
              weights=weights, followup=followup, entry=entry, exit=exit,
              obswin=obswin, scheme=scheme,
              model=model, prior.var=scale^2, post.var=post.var,method=method,
              mtd=rec, include=include, pid=pid, model.detail=model.detail,intcpt=intcpt,
              ptox=ptox, ptoxL=ptoxL, ptoxU=ptoxU, conf.level=conf.level,
              patient.detail=patient.detail,tite=tite,dosescaled=dosescaled)
  class(foo) <- "mtd"
  foo
}
