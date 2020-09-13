titesim <- function(PI, # true toxicity probabilites
                    prior, # initial guess of toxicity probabilites
                    target, # target DLT rate
                    n, # number of patients
                    x0, # starting dose level for 1-stage CRM
                    nsim=1, # number of simulations
                    restrict=TRUE, # safety !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    obswin=1,
                    tgrp=obswin,
                    rate=1,
                    accrual="fixed",
                    surv="uniform",
                    scheme="linear",
                    count=TRUE,
                    method="bayes",
                    model="empiric",
                    intcpt=3,
                    scale=sqrt(1.34),
                    seed=1009) {
  if (nsim==1) {
    foo <- onetite(PI,prior,target,n,x0,restrict=restrict,obswin=obswin,tgrp=tgrp,rate=rate,accrual=accrual,
                   surv=surv,scheme=scheme,method=method,model=model,intcpt=intcpt,scale=scale,seed=seed)
  }
  else {
    foo <- mtite(nsim,PI,prior,target,n,x0,restrict=restrict,obswin=obswin,tgrp=tgrp,rate=rate,accrual=accrual,count=count,
                 surv=surv,scheme=scheme,method=method,model=model,intcpt=intcpt,scale=scale,seed=seed)
  }
  foo
}



