fit.mix.model <- function(means, markers, tau, n.states = 3, common.sigma2 = FALSE, tol = 1e-6) {
  n <- length(means)
  temp <- apply(matrix(runif(n * (n.states - 1)), nrow = n.states - 1, ncol = n), 2, sort)
  p.states <- t(rbind(temp, 1) - rbind(0, temp))
  #init.states <- unclass(cut(rank(means), breaks = n.states))
  #p.states <- matrix(0, nrow = n, ncol = n.states)
  #for (i in 1:n) p.states[i, init.states[i]] <- 1
  change <- Inf
  loglik <- -Inf
  taui <- tau/markers
  sigma2 <- rep(1, n.states)
  while (change > tol) {
    ## M-step
    p0 <- colSums(p.states)/sum(p.states)
    mu <- colSums(means * p.states/(outer(taui, sigma2, "+")))/colSums(p.states/(outer(taui, sigma2, "+")))
    target.each <- list()
    gen.tar <- function(i) {
      f <- function(s2) {
        return(sum((means - mu[i])^2 * p.states[, i]/(s2 + taui)^2) - sum(p.states[, i]/(s2 + taui)))
      }
      return(f)
    }
    for(i in 1:n.states) {
      target.each[[i]] <- gen.tar(i)
    }
    if(common.sigma2) {
      target.common <- function(s2) {
        res <- 0
        for(i in 1:n.states) {
          res <- res + target.each[[i]](s2)
        }
        return(res)
      }
      if (target.common(0) < 0) {
        sigma2 <- rep(0, n.states)
      } else {
        sigma2 <- rep(uniroot(target.common, c(0, 100))$root, n.states)
        #sigma2 <- rep(sum(outer(means, mu, "-")^2 * markers * p.states)/sum(p.states), n.states)
      }
    } else {
      for(i in 1:n.states) {
        if (target.each[[i]](0) < 0) {
          sigma2[i] <- 0
        } else {
          sigma2[i] <- uniroot(target.each[[i]], c(0, 100))$root
        }
      }
      #sigma2 <- colSums(outer(means, mu, "-")^2 * markers * p.states)/colSums(p.states)
    }
    ## E-step
    logprob <- matrix(dnorm(rep(means, n.states), rep(mu, each = n), sqrt(outer(taui, sigma2, "+")), log = TRUE), ncol = n.states) + rep(log(p0), each = n)
    #logprob <- matrix(dnorm(rep(means, n.states), rep(mu, each = n), sqrt((1/markers) %o% sigma2), log = TRUE), ncol = n.states) + rep(log(p0), each = n)
    norm.fac <- apply(logprob, 1, max)
    prob1 <- exp(logprob - norm.fac)
    p.states <- prob1/apply(prob1, 1, sum)
    ## Log-likelihood
    old.loglik <- loglik
    loglik <- sum(log(rowSums(prob1)) + norm.fac)
    change <- loglik - old.loglik
    ##print(change)
  }
  return(list(loglik = loglik, p0 = p0, mu = mu, sigma2 = sigma2, p.states = p.states))
}
