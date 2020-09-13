crmh = function(a,x,y,w,s) {  ## posterior
  v = exp(-a^2/2/s^2)
  for (i in 1:length(x)) {
    v = v * ((x[i]^exp(a))^y[i])*(((1-w[i]*x[i]^exp(a))^(1-y[i])))
  }
  return(v)
}

crmht = function(a,x,y,w,s) { ## posterior times x
  v = a * exp(-a^2/2/s^2)
  for (i in 1:length(x)) {
    v = v * ((x[i]^exp(a))^y[i])*(((1-w[i]*x[i]^exp(a))^(1-y[i])))
  }
  return(v)
}

crmht2 = function(a,x,y,w,s) { ## posterior times x^2
  v = a^2 * exp(-a^2/2/s^2)
  for (i in 1:length(x)) {
    v = v * ((x[i]^exp(a))^y[i])*(((1-w[i]*x[i]^exp(a))^(1-y[i])))
  }
  return(v)
}


crmhlgt <- function(a,x,y,w,s,alp0)  { ## posterior logit model
  v = exp(-a^2/2/s^2)
  for (i in 1:length(x)) {
    PSI <- (1 + exp(-alp0-exp(a)*x[i]))^{-1}
    v <- v * (PSI^y[i]) * (1-w[i]*PSI)^(1-y[i])
  }
  return(v)
}
crmhtlgt <- function(a,x,y,w,s,alp0)  { ## posterior times x
  v = a * exp(-a^2/2/s^2)
  for (i in 1:length(x)) {
    PSI <- (1 + exp(-alp0-exp(a)*x[i]))^{-1}
    v <- v * (PSI^y[i]) * (1-w[i]*PSI)^(1-y[i])
  }
  return(v)
}
crmht2lgt <- function(a,x,y,w,s,alp0)  { ## posterior times x^2
  v = a^2 * exp(-a^2/2/s^2)
  for (i in 1:length(x)) {
    PSI <- (1 + exp(-alp0-exp(a)*x[i]))^{-1}
    v <- v * (PSI^y[i]) * (1-w[i]*PSI)^(1-y[i])
  }
  return(v)
}



lcrm <- function(a,x,y,w) { #loglikelihood of empiric function
  v <- 0
  for (i in 1:length(x))
    v <- v + y[i]*log(x[i])*exp(a) + (1-y[i])*log(1 - w[i]*x[i]^exp(a))
  return(v)
}
lcrmlgt <- function(a,x,y,w,alp0) { #loglikelihood of logit function
  v <- 0
  for (i in 1:length(x)) {
    PSI <- (1 + exp(-alp0-exp(a)*x[i]))^{-1}
    v <- v + y[i]*log(PSI) + (1-y[i])*log(1-w[i]*PSI)
  }
  return(v)
}
