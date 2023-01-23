# Exact test --------------------------------------------------------------
r.freq= function(theta0, n, alpha = 0.05, hp)
{
  ris.max =  qpois(1 - 10^-8, n*theta0)
  ris = 0:ris.max 
  if(hp == "less")
  {
    prob.ris.cum = ppois(ris,n*theta0)  # probability to observe a smaller value
    r = max(ris[prob.ris.cum<=alpha])# Critical value 
  }
  else
  {
    prob.ris = dpois(rev(ris),n*theta0)  # Probability mass for each result
    prob.ris.cum = cumsum(prob.ris)  # probability to observe a larger result
    r = min(rev(ris)[prob.ris.cum<=alpha])# Critical value 
  }
  return(r)
}
r.freq.v = Vectorize(r.freq, vectorize.args = "n")



r.bayes = function(theta0, n, alphaA, betaA, lambda = 0.95, hp)
{
  ris.max =  qpois(1 - 10^-8, n*theta0)
  ris = 0:ris.max
  al = alphaA + ris
  be = betaA + n
  
  if(hp == "less")
  {
    bayes.cond = pgamma(theta0, shape = al, rate = be)
    r.bayes = max(ris[bayes.cond>lambda])
    return(r.bayes)
  }
  else{
    bayes.cond = 1 - pgamma(theta0, shape = rev(al), rate = be)
    r.bayes = min(rev(ris)[bayes.cond>lambda])
    return(r.bayes)
  }
}
r.bayes.v = Vectorize(r.bayes, vectorize.args = "n")

# Conditional and predictive approach -------------------------------------

## Conditional power
power.c = function(thetaD, n, r, hp)
{
  if(hp == "less"){
    power = ppois(r,n*thetaD)
  }
  else
  {
    power = 1 - ppois(r-1, n*thetaD)
  }
  return(power)
}

## Predictive power
power.p = function(alphaD, betaD, n, r, hp)
{
  if(hp == "less"){
    power = pnbinom(r,alphaD, betaD/(betaD + n))
  }
  else{
    power = 1 - pnbinom(r-1,alphaD, betaD/(betaD + n))
  }
  return(power)
}



# Power functions ---------------------------------------------------------
eta = function(theta0, hp, thetaD, 
               alphaD, betaD, 
               alphaA, betaA, 
               n, alpha, lambda, 
               power.type)
{
  
  r <- 0
  power <- 0
  
  ## Critical value calculation
  if(grepl("F",power.type)){
    r = r.freq(theta0 = theta0, n = n, alpha = alpha, hp = hp)
  }
  else{
    r = r.bayes(theta0 = theta0, n = n, alphaA = alphaA, betaA = betaA, lambda = lambda, hp = hp)
  }
  ## Approach to eliminate dependency form theta
  if(grepl("C", power.type))
  {
    power <- power.c(thetaD = thetaD, n = n, r = r, hp = hp)
  }
  else{
    power <- power.p(alphaD = alphaD, betaD = betaD, n = n, r = r, hp = hp)
  }
  
  return(list(power, r))
}
eta.v = Vectorize(eta, vectorize.args = "n")



# SSD determination -------------------------------------------------------
n_min_opt = function(dim.n, vett.prob, r.vett, lev.pot = 0.8)
{
  # Classic criterion
  n.cl = min(dim.n[vett.prob>lev.pot])
  r.cl = r.vett[dim.n == n.cl]
  # Conservative criterion
  n.co = r.co = NA
  if (min(vett.prob)<=lev.pot && max(vett.prob)>=lev.pot){
    n.co = dim.n[max(which(vett.prob<=lev.pot))]+1 
    r.co = r.vett[dim.n == n.co]
  }
  if (min(vett.prob)>lev.pot) { 
    n.co = 0
    r.co = 0  }
  if (max(vett.prob)<lev.pot) { 
    n.co = NA
    r.co = NA  }
  return(list(n.cl, r.cl, n.co, r.co)) 
}


# Prior Parameters --------------------------------------------------------
# Link between prior mode and sample size and alpha and beta

prior.pois = function(theta0, n0)
{
  alpha = n0*theta0 + 1
  beta = n0
  return(c(alpha = alpha,beta = beta))
}

diff.int.pois = function(n0, thetaM, l.inf, l.sup, lev = 0.95)
{
  priors = prior.pois(thetaM,n0)
  prob = integrate(dgamma, l.inf, l.sup, shape = priors[1], rate = priors[2])
  diff = pgamma(l.sup, shape = priors[1], rate = priors[2]) - pgamma(l.inf, shape = priors[1], rate = priors[2]) - lev
  return(diff)
}


n0.star.pois = function(thetaM, l.inf, l.sup, lower = 1, upper = 10000, lev = 0.95)
{
  n0 = uniroot(diff.int.pois, lower = lower, upper = upper, thetaM = thetaM, l.inf = l.inf, l.sup = l.sup, lev = lev)$root
  return(n0)
}

hyper.gamma.int = function(thetaM, lev = 0.95, l.inf, l.sup, lower = 0.05, upper = 1000){
  n0 = n0.star.pois(thetaM, l.inf, l.sup, lower, upper, lev)
  return(c(prior.pois(thetaM, n0), n0 = n0))
}

