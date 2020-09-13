mtite <- function(nsim, PI, prior, target, n, x0, restrict=TRUE, obswin=1, tgrp=obswin, rate=1, accrual="fixed",
                  surv="uniform", scheme="linear", count=TRUE, method="bayes",model="empiric",intcpt=3,scale=sqrt(1.34),seed=1009) {
  set.seed(seed)
  if (nsim<=1) { stop(" nsim less than 2!  For single simulation replicate, use function `onetite'."); }
  nexpt <- ntox <- sel <- rep(0,length(prior))
  BETAHAT <- matrix(rep(NA,nsim*n), nrow=nsim)
  final.est <- DURATION <- rep(NA,nsim)

  for (r in 1:nsim) {
    if (count) { cat("simulation number:",r,"\n"); }

    if (accrual=="fixed") { next.arrival <- obswin/rate; }
    else if (accrual=="poisson") { next.arrival <- rexp(1, rate/obswin); }

    if (length(x0)>1) {
      if (length(x0)!=n) { stop(" Initial design has a different sample size than that is specified!"); }
      bethat <- rep(0,n)
      u <- y <- level <- arrival <- rep(NA,n)
      m <- 1
      while (TRUE) {
        level[m] <- cur <- x0[m]
        if (is.na(arrival[m])) { arrival[m] <- next.arrival; }
        if (is.na(y[m])) {
          if (surv=="uniform") {
            y[m] <- ynew <- rbinom(1,1,PI[cur])
            if (ynew) { unew <- runif(1,0,obswin); }
            else { unew <- Inf; }
            u[m] <- unew
            utox <- u + arrival
          }
        }
        if (accrual=="fixed") { next.arrival <- next.arrival + obswin/rate; }
        else if (accrual=="poisson") { next.arrival <- next.arrival + rexp(1, rate/obswin); }
        B <- rep(0,n); B[utox<=next.arrival] <- 1;
        if (sum(B)>0 | m==(n-1)) { break; }
        if (x0[m+1]==cur | (next.arrival-arrival[m])>=tgrp) { m <- m+1; }
      }
      if (m==(n-1)) {
        if (sum(B)==0) { cur <- x0[n]; }
        else {
          censor <- pmin(next.arrival, utox) - arrival;
          followup <- pmin(censor, obswin)
          if (mean(B[1:m])==1) {
            obj <- titecrm(prior, target, B[1:m], level[1:m], followup=followup[1:m],obswin=obswin,scheme=scheme,model=model,
                           intcpt=intcpt,scale=500,var.est=FALSE)
          }
          else {
            obj <- titecrm(prior, target, B[1:m], level[1:m], followup=followup[1:m],obswin=obswin, scheme=scheme,method=method,model=model,
                           intcpt=intcpt,scale=scale,var.est=FALSE)
          }
          if (restrict) { cur <- min(obj$mtd, (cur+1)); }
          else { cur <- obj$mtd; }
          bethat[m+1] <- obj$est
        }
        arrival[n] <- next.arrival
        level[n] <- cur
        if (surv=="uniform") {
          y[n] <- ynew <- rbinom(1,1, PI[cur])
          if (ynew) { unew <- runif(1,0,obswin); }
          else { unew <- Inf; }
          u[n] <- unew; utox <- u + arrival
        }
      }
      else {
        censor <- pmin(next.arrival, utox) - arrival
        followup <- pmin(censor, obswin)
        if (mean(B[1:m])==1) {
          obj <- titecrm(prior, target, B[1:m], level[1:m], followup=followup[1:m],obswin=obswin,scheme=scheme,model=model,
                         intcpt=intcpt,scale=500,var.est=FALSE)
        }
        else {
          obj <- titecrm(prior, target, B[1:m], level[1:m], followup=followup[1:m],obswin=obswin, scheme=scheme,method=method,model=model,
                         intcpt=intcpt,scale=scale,var.est=FALSE)
        }
        if (restrict) { cur <- min(obj$mtd, (cur+1)); }
        else { cur <- obj$mtd; }
        bethat[m+1] <- obj$est
        #        if (accrual=="fixed") { next.arrival <- next.arrival + obswin/rate; }
        #        else if (accrual=="poisson") { next.arrival <- next.arrival + rexp(1, rate/obswin); }

        for (i in (m+1):(n-1)) {
          arrival[i] <- next.arrival; level[i] <- cur;
          if (surv=="uniform") {
            y[i] <- ynew <- rbinom(1,1,PI[cur])
            if (ynew) unew <- runif(1,0,obswin)
            else unew <- Inf
            u[i] <- unew; utox <- u + arrival;
          }
          if (accrual=="fixed") next.arrival <- next.arrival + obswin/rate
          else if (accrual=="poisson") next.arrival <- next.arrival + rexp(1,rate/obswin)
          B <- rep(0,n);  B[utox<=next.arrival] <- 1;
          censor <- pmin(next.arrival,utox) - arrival; followup <- pmin(censor,obswin);
          if (mean(B[1:i])==1) {
            obj <- titecrm(prior, target, B[1:i], level[1:i], followup=followup[1:i],obswin=obswin,scheme=scheme,model=model,
                           intcpt=intcpt,scale=500,var.est=FALSE)
          }
          else {
            obj <- titecrm(prior, target, B[1:i], level[1:i], followup=followup[1:i],obswin=obswin, scheme=scheme,method=method,model=model,
                           intcpt=intcpt,scale=scale,var.est=FALSE)
          }
          if (restrict) { cur <- min(obj$mtd, (cur+1)); }
          else { cur <- obj$mtd; }
          bethat[i+1] <- obj$est
        }
        BETAHAT[r,] <- bethat

        arrival[n] <- next.arrival; level[n] <- cur;
        if (surv=="uniform") {
          y[n] <- ynew <- rbinom(1,1,PI[cur])
          if (ynew) unew <- runif(1,0,obswin)
          else unew <- Inf
          u[n] <- unew; utox <- u + arrival
        }
      }
      if (method=="mle" & (sum(y)==0|mean(y)==1)) {
        finalobj <- titecrm(prior,target,y,level,weights=rep(1,n),model=model,intcpt=intcpt,scale=500,var.est=FALSE)
        msg <- "Warning: mle is approximated"
      }
      else {
        finalobj <- titecrm(prior, target, y, level, weights=rep(1,n),method=method,model=model,intcpt=intcpt,scale=scale,var.est=FALSE)
        msg <- "Okay"
      }
      cur <- finalobj$mtd
      est <- finalobj$est
    }
    else {
      if (method=="mle") { stop(" Require an initial design for mle-CRM!"); }
      bethat <- 0
      u <- y <- level <- arrival <- NULL
      cur <- x0

      for (i in 1:(n-1)) {
        arrival <- c(arrival, next.arrival)
        level <- c(level, cur)
        if (surv=="uniform") {
          ynew <- rbinom(1,1,PI[cur])
          if (ynew) unew <- runif(1,0,obswin)
          else unew <- Inf
          y <- c(y,ynew); u <- c(u,unew); utox <- u + arrival;
        }
        if (accrual=="fixed") { next.arrival <- next.arrival + obswin/rate; }
        else if (accrual=="poisson") { next.arrival <- next.arrival + rexp(1, rate/obswin); }
        B <- rep(0,length(y));  B[utox<=next.arrival] <- 1;
        censor <- pmin(next.arrival,utox) - arrival
        followup <- pmin(censor,obswin)
        obj <- titecrm(prior,target,B,level,followup=followup,obswin=obswin, scheme=scheme, method=method,model=model,intcpt=intcpt,scale=scale,var.est=FALSE)
        if (restrict) {
          cur <- min(obj$mtd, (cur+1))
        }
        else {
          cur <- obj$mtd
        }
        bethat <- c(bethat, obj$est)
      }
      BETAHAT[r,] <- bethat
      arrival <- c(arrival, next.arrival)
      level <- c(level,cur)
      if (surv=="uniform") {
        ynew <- rbinom(1,1,PI[cur])
        if (ynew) unew <- runif(1,0,obswin)
        else unew <- Inf
        y <- c(y,ynew); u <- c(u,unew); utox <- u + arrival;
      }
      finalobj <- titecrm(prior,target,y,level,weights=rep(1,n), method=method,model=model,intcpt=intcpt,scale=scale,var.est=FALSE)
      cur <- finalobj$mtd
      est <- finalobj$est
      msg <- "Okay"
    }

    sel[cur] <- sel[cur] + 1
    final.est[r] <- est
    DURATION[r] <- max(arrival) + obswin
    for (k in 1:length(prior)) {
      nexpt[k] <- nexpt[k] + length(which(level==k))
      ntox[k] <- ntox[k] + length(which(y==1 & level==k))
    }
  }
  sel <- sel/nsim
  nexpt <- nexpt/nsim
  ntox <- ntox/nsim
  if (length(x0)==1) { design <- paste("TITE-CRM starting at dose",x0); }
  else { design <- "Two-stage TITE-CRM"; }

  foo <- list(PI=PI, prior=prior, target=target, n=n, x0=x0, nsim=nsim,
              MTD=sel, level=nexpt, tox=ntox, beta.hat=BETAHAT, final.est=final.est,
              Duration=DURATION,
              design=design, method=method, prior.var=scale^2, model=model, intcpt=intcpt,
              restriction=restrict, seed=seed, tite=TRUE, dosescaled=finalobj$dosescaled, msg=msg,
              obswin=obswin, tgrp=tgrp, rate=rate, accrual=accrual, scheme=scheme)

  class(foo) <- "sim"
  foo
}
