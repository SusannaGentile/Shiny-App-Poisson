# Section 2 ---------------------------------------------------------------

# We list the main inputs, which are common to the majority of the functions.
# The names of the parameters are as coherent as possible with the notation
# adopted through the article.

# Input descriptions:
# - theta0: event rate under the null hypothesis
# - hp: alternative hypothesis to be tested, "less" or "greater"
# - n: sample size
# - k: critical value, regardless of the analysis framework
# - alpha: first error rate
# - epsilon: we define the Bayesian significance threshold as 1-epsilon
# - alphaA, betaA: analysis prior hyperparameters
# - thetaD: design value
# - alphaD, betaD: design prior hyperparameters
# - gamma: desired power level
# - power.type: a two-characters string indicating the analysis framework
#               ("F" or "B") and the power computation method ("C" or "P")
# Some functions have specific inputs, which are listed in their description

##### Critical value computation

# Function to compute the frequentist critical value (Section 2.1)
# Output: critical value of the exact frequentist test

k.freq = function(theta0, n, alpha = 0.05, hp = "less") {
  ris.max = qpois(1 - 10^(-10), n * theta0)
  ris = switch(hp,
    "greater" = ris.max:0,
    "less" = 0:ris.max
  )
  prob.ris = dpois(ris, n * theta0)
  prob.ris.cum = cumsum(prob.ris)
  k.freq = switch(hp,
    "greater" = min(ris[prob.ris.cum <= alpha]),
    "less" = max(ris[prob.ris.cum <= alpha])
  )
  return(k.freq)
}

k.freq.v = Vectorize(k.freq, vectorize.args = "n")

# Function to compute the Bayesian critical value (Section 2.3)
# Output: critical value in a Bayesian perspective

k.bayes = function(theta0, n, hp = "less",
                    alphaA, betaA, epsilon = 0.05) {
  ris.max = qpois(1 - 10^(-10), n * theta0)
  ris = switch(hp,
    "greater" = ris.max:0,
    "less" = 0:ris.max
  )
  al = alphaA + ris
  be = betaA + n
  pgamma.theta0 = pgamma(theta0, shape = al, rate = be)
  k.bayes = switch(hp,
    "greater" = min(ris[(1 - pgamma.theta0) > 1 - epsilon]),
    "less" = max(ris[pgamma.theta0 > 1 - epsilon])
  )
    return(k.bayes)
}

k.bayes.v = Vectorize(k.bayes, vectorize.args = "n")

##### Conditional and predictive approach

# Function implementing the conditional approach to assume that theta belongs
# to the alternative hypothesis
# Output: probability of rejecting the null hypothesis with respect to the
# sampling distribution of the test statistic Sn = pois(n*thetaD)

conditional = function(thetaD, n, k, hp = "less") {
  if (hp == "less") {
    power = ppois(k, n * thetaD)
  } else {
    power = 1 - ppois(k - 1, n * thetaD)
  }
  return(power)
}

# Function implementing the predictive approach to assume that theta belongs
# to the alternative hypothesis
# Output: probability of rejecting the null hypothesis with respect to the
# predictive distribution of Sn = bin-neg(n, alphaD, betaD/(betaD + n)).

predictive = function(alphaD, betaD, n, k, hp = "less") {
  if (hp == "less") {
    power = pnbinom(k, alphaD, betaD / (betaD + n))
  } else {
    power = 1 - pnbinom(k - 1, alphaD, betaD / (betaD + n))
  }
  return(power)
}

##### Power computation

# The function consists of two step depending on the "power.type" argument.
# First, it computes the critical value according to the first letter, then
# the probability of rejecting H_0 based on the second letter.

# Output:
# - the value of the required power function
# - if the output.type input is modified, the critical value, either
#   frequentist or bayesian depending on the first letter of the
#   power.type argument

eta = function(n, theta0, hp = "less", power.type,
                thetaD = NULL, alphaD = NULL, betaD = NULL,
                alphaA = NULL, betaA = NULL,
                alpha = 0.05, epsilon = 0.05,
                output.type = "power") {
  if (!(power.type %in% c("FC", "FP", "BC", "BP"))) {
    stop("The power type should be among FC, FP, BC, BP")
  }

  k = 0
  power = 0

  ## Critical value calculation
  if (grepl("F", power.type)) {
    k = k.freq(theta0 = theta0, n = n, alpha = alpha, hp = hp)
  } else {
    k = k.bayes(
      theta0 = theta0, n = n, alphaA = alphaA, betaA = betaA,
      epsilon = epsilon, hp = hp
    )
  }

  ## Approach to eliminate dependency from theta
  if (grepl("C", power.type)) {
    power = conditional(
      thetaD = thetaD,
      n = n, k = k, hp = hp
    )
  } else {
    power = predictive(
      alphaD = alphaD, betaD = betaD,
      n = n, k = k, hp = hp
    )
  }
  ## Which output to return

  if (output.type == "power") {
    return(power)
  } else {
    return(c(power, k))
  }
}

# Vectorization with respect to the sample size
eta.v = Vectorize(eta, vectorize.args = "n")

###### SSD power-based criteria

# Function which computes the optimal sample size at level gamma according to:
# - the classic criterion: the optimal n is the minimum n so that the power
#   exceeds gamma
# - the conservative criterion: it accounts for the saw-toothed behavior by
#   selecting as optimal sample size the minimum n so that the power exceeds
#   gamma for every sample size greater than or equal to the optimal one

# The inputs are:
# - dim.n: a vector of candidates for the optimal sample size
# - power.vec: a vector a probabilities representing the power associated to
#   every sample size
# - k.vec: a vector containing the critical values associated to
#          every sample size
# The implementation of the two criteria does not depend on which power function
# we consider. Thus, the power.vec and k.vec are still generic

# Output: a list of four elements
# [[1]] sample size according to the classic criterion
# [[2]] critical value when the sample size is the optimal one according to the
#       classic criterion
# [[3]] sample size according to the conservative criterion
# [[4]] critical value when the sample size is the optimal one according to the
#       conservative criterion

n.min.opt = function(dim.n, power.vec, k.vec, gamma = 0.8) {
  # Classic criterion
  n.cl = min(dim.n[power.vec > gamma])
  k.cl = k.vec[dim.n == n.cl]

  # Conservative criterion
  n.co = r.co = NA
  if (min(power.vec) <= gamma && max(power.vec) >= gamma) {
    n.co = dim.n[max(which(power.vec <= gamma))] + 1
    k.co = k.vec[dim.n == n.co]
  }

  if (min(power.vec) > gamma) {
    n.co = 0
    k.co = 0
  } # The sample size is fulfilled for every sample size
  if (max(power.vec) < gamma) {
    warning("nmax is too small to reach the desired power level")
    n.co = NA
    k.co = NA
  } #
  return(list(n.cl, k.cl, n.co, k.co))
}

# Function combining the eta and the n.min.opt functions.
# It computes the inputs of the n.min.opt function on using the eta function,
# then it calls the n.min.opt function
# Specific inputs:
# - n.min, n.max, step: arguments to define the dim.n input of the n.min.opt
#                       function. The step argument allows to consider also not 
#                       integer sample sizes
# - nA, thetaA: prior sample size and prior mode of the analysis prior
# - nD: prior sample size for the design prior

# Output: a list of 7 elements
# [[1]] sample size according to the classic criterion
# [[2]] critical value when the sample size is the optimal one according to the
#       classic criterion
# [[3]] sample size according to the conservative criterion
# [[4]] critical value when the sample size is the optimal one according to the
#       conservative criterion
# [[5]] the vector of candidates for the optimal sample size
# [[7]] a vector containing the critical values associated to
#       every sample size
# [[6]] a vector containing the power values associated to every sample size

n.pois = function(theta0, hp = "less", power.type = "FC", gamma = 0.8,
                   thetaD = NULL, nD = NULL,
                   alphaD = thetaD * nD + 1, betaD = nD,
                   thetaA = NULL, nA = NULL,
                   alphaA = thetaA * nA + 1, betaA = nA,
                   alpha = 0.05, epsilon = 0.05,
                   n.max = 300, n.min = 5, step = 1) {
  dim.n = seq(from = n.min, to = n.max, by = step)
  power = eta.v(
    theta0 = theta0, hp = hp, thetaD = thetaD, alphaD = alphaD,
    betaD = betaD, alphaA = alphaA, betaA = betaA, n = dim.n,
    alpha = alpha, epsilon = epsilon, power.type = power.type,
    output.type = "complete"
  )

  power.vec = power[1, ]
  k.vec = power[2, ]
  results = n.min.opt(
    dim.n = dim.n, power.vec = power.vec, k.vec = k.vec,
    gamma = gamma
  )

  return(append(results, list(dim.n, k.vec, power.vec)))
}

# We vectorize the function
n.pois.vec = Vectorize(n.pois)


# Section 3 ------------------------------------------------
# These functions implement two possible strategies to elicit the priors
# distributions.

# Inputs:
# thetaM: prior mode
# n0: prior sample size
# l.inf, l.sup: interval to which the prior should assign a specified probability
# lev: probability to assign to the interval [l.inf, l.sup]
# lower, upper: see uniroot

# Hyperparameters given the prior mode and prior sample size
# Output: hyperparameters alpha and beta
hyperparameters = function(thetaM, n0) {
  alpha = n0 * thetaM + 1
  beta = n0
  return(c(alpha = alpha, beta = beta))
}

hyper.vec = Vectorize(hyperparameters)

# Function computing the difference between the probability assigned to
# [l.inf, l.sup] by the prior distribution with prior mode thetaM and prior
# sample size n0 and the target probability lev
diff.int = function(n0, thetaM, l.inf, l.sup, lev = 0.999) {
  prior = hyperparameters(thetaM, n0)
  prob = integrate(dgamma, l.inf, l.sup, shape = prior[1], rate = prior[2])
  diff = prob$value - lev
  return(diff)
}

# Function which searches the interval from lower to upper for a root of the
# function diff.int with respect to the prior sample size n0
# Output: optimal n0
n0.star = function(thetaM, l.inf, l.sup, lev = 0.95,
                    lower = 1, upper = 1000) {
  n0 = uniroot(diff.int,
    lower = lower, upper = upper, thetaM = thetaM,
    l.inf = l.inf, l.sup = l.sup, lev = lev
  )$root
  return(n0)
}

n0.star.vec = Vectorize(n0.star)

# Hyperparameters given the prior mode, the interval [l.inf, l.sup] and the
# probability lev
# Output: hyperparameters alpha and beta
hyper.int = function(thetaM, lev = 0.95, l.inf, l.sup, lower = 1,
                      upper = 1000) {
  n0 = n0.star(
    thetaM = thetaM, lev = lev,
    l.inf = l.inf, l.sup = l.sup,
    lower = lower, upper = upper
  )
  prior = hyperparameters(thetaM = thetaM, n0 = n0)
  return(c(prior, n0 = n0))
}

hyper.int.vec = Vectorize(hyper.int)