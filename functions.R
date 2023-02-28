##' Count the number of observations for the parent 
##' 
##' @param child_count vector of real, indicating the number of 
##' observations falling into the child
##' 
##' @return vector of length 1/2*length(child_count), the 
##' number of observations falling into the parents
get_parent_count <- function(child_count){
  parent_len <- length(child_count) / 2
  return (child_count[2*(1:parent_len)-1] + child_count[2*(1:parent_len)])
}


##' Calculate the updated/posterior 
##' probability measure for the partitions
##' 
##' @param theta vector of parameter theta
##' @param maxK scalar, maximum tree levels i.e.depth of the tree
##' 
##' @return vector of posterior probabilities for each partition
Cal_prob <- function(theta, maxK){
  prob <- c(theta[1], 1-theta[1])
  prob_i <- c(theta[1], 1-theta[1])
  for(i in 2:maxK){
    new_prob_i <- c()
    for(j in 1:length(prob_i)){
      x <- theta[2^(i-1) + j - 1]
      new_prob_i <- c(new_prob_i, prob_i[j]*x, prob_i[j]*(1-x))
    }
    prob_i <- new_prob_i
    prob <- c(prob, new_prob_i)
  }
  return (prob)
}

##' Calculate the log marginal probability of 
##' a set of observations y
##' 
##' @param alpha vector of parameter alpha
##' @param maxK scalar, maximum tree levels 
##' 
##' @return scalar, the log marginal probability of y
log_Marginal_prob <- function(a_prior,a_post){
  log_marg <- 0
  for(i in 1:(length(a_prior)/2)){
    log_marg <- log_marg + (lgamma(a_post[2*i-1])+lgamma(a_post[2*i])+
      lgamma(a_prior[2*i-1]+a_prior[2*i])-
      (lgamma(a_post[2*i-1]+a_post[2*i])+
         lgamma(a_prior[2*i-1])+lgamma(a_prior[2*i])))
  }
  return (log_marg)
}


##' Draw theta samples
##' 
##' @param alpha vector of alpha
##' 
##' @return vector of theta drawing from a specified beta
draw_theta <- function(alpha=a_post){
  theta <- c()
  for(i in 1:(length(alpha)/2)){
    theta <- c(theta, rbeta(1, alpha[2*i-1], alpha[2*i]))
  }
  return (theta)
}


##' Perform a Polya Tree model 
##' 
##' @param a scalar, left point of the bounded support
##' @param b scalar, right point of the bounded support
##' @param maxK scalar, maximum level of the tree (i.e.depth)
##' @param a_prior vector, prior distribution for alpha
##' @param y_obs vector of observations
##' 
##' @return vector, the posterior probability measure after observing data
PT_update <- function(a, b, maxK=10, a_prior, y_obs){
  # Total number of nodes in the tree
  nParam <- 2^(maxK + 1) - 2
  
  # Break the whole set into disjoint partitions (leaf)
  breaks <- a + (1:(2^maxK-1)) / (2^maxK) * (b-a) 
  
  # Make all observations fall into the bounded support
  y_obs[y_obs < a] <- a
  y_obs[y_obs > b] <- b
  
  # Count the number of observations falling into 
  # the partitions (nodes in the bottom level of the tree)
  leaf_count <- hist(y_obs, breaks=c(-Inf, breaks, Inf), plot=FALSE)$counts
  
  # Find the number of observations falling into each nodes recursively
  count <- c()
  child_count <- leaf_count
  while(length(child_count) > 1){
    count <- c(child_count, count)
    parent_count <- get_parent_count(child_count)
    child_count <- parent_count
  }
  
  # Update alpha posterior
  a_post <- a_prior + count
  
  # Draw samples from the posterior of theta
  theta <- draw_theta(a_post)
  
  prob <- Cal_prob(theta, maxK)
  
  return (list(prob=prob,a_post=a_post))
}

##' Generate alpha priors which results absolute continuous
##' random probability measures -- alpha=cm^2
##' 
##' @param c scalar, hyperparameter 
##' @param maxK scalar, maximum tree levels 
##' @param G0 centering distribution
##' i.e.depth of the tree
##' 
##' @return vector of alpha prior
alpha_priors <- function(c=1, maxK){
  alpha <- c()
  for(m in 1:(2^maxK-1)){
    alpha[2*m-1] <- c*m^2
    alpha[2*m] <- c*m^2
  }
  return (alpha)
}

##' Find change point location
##' 
##' @param y_obs vector of ordered data
##' @param maxK scalar, maximum tree levels 
##' 
##' @return scalar, the changepoint location
find_change_location <- function(y_obs, maxK=7, a_prior){
  log_marg_prev <- rep(-Inf, len=length(y_obs))
  log_marg_after <- rep(-Inf, len=length(y_obs))
  n_obs <- length(y_obs)
  ign <- n_obs/50
  for (tau in ign:(n_obs-ign)){
    y_obs_prev <- y_obs[1:tau-1]
    y_obs_after <- y_obs[tau:n_obs]
    
    a1 <- floor(min(y_obs_prev))
    b1 <- ceiling(max(y_obs_prev))
    a_post_prev <- PT_update(a1,b1,maxK,a_prior, y_obs_prev)$a_post
    
    a2 <- floor(min(y_obs_after))
    b2 <- ceiling(max(y_obs_after))
    a_post_after <- PT_update(a2,b2,maxK,a_prior, y_obs_after)$a_post
    
    log_marg_prev[tau] <- log_Marginal_prob(a_prior,a_post_prev)
    log_marg_after[tau] <- log_Marginal_prob(a_prior,a_post_after)
  }
  log_marg <- log_marg_prev + log_marg_after
  change <- which.max(log_marg)
  return(change)
}

##' Compute the covering metric 
##' 
##' @param true_idx scalar, the true change location
##' @param est_idx scalar, the estimated change location
##' @param vec_len the length of the underlying time series
##' 
##' @return scalar, the changepoint location
covering_metric <- function(true_idx, est_idx, vec_len){
  d <- abs(true_idx-est_idx)
  l <- min(true_idx, vec_len-true_idx)
  return(d/l)
}

##' Nonparametric Bayesian hypothesis test
##' Calculate the Bayes factor for H0 against H1
##' 
##' @param log_marg.H0 scalar, the log marginal probability under H0
##' @param log_marg.H1 vector, log marginal under H1 across 
##' every changepoint configuration
##' 
##' @return scalar, Bayes factor for the hypothesis test
cal_BF <- function(a_prior, a_post.H0, a_post_prev,a_post_after){
  b <- c()
  for(i in 1:(length(a_prior)/2)){
    b[i] <- lgamma(a_prior[2*i-1])+lgamma(a_prior[2*i])+
      lgamma(a_post.H0[2*i-1])+lgamma(a_post.H0[2*i])-
      lgamma(a_prior[2*i-1]+a_prior[2*i])-
      lgamma(a_post.H0[2*i-1]+a_post.H0[2*i])+
      lgamma(a_post_prev[2*i-1]+a_post_prev[2*i])+
      lgamma(a_post_after[2*i-1]+a_post_after[2*i])-
      lgamma(a_post_prev[2*i-1])-
      lgamma(a_post_prev[2*i])-
      lgamma(a_post_after[2*i-1])-
      lgamma(a_post_after[2*i])
    
  }
  log_odds_ratio <- sum(b)
  odds.post.H0 <- 1/(1+exp(-log_odds_ratio))
  return(sum(b))
}



##' Exponential kernel 
##' @param X1 numeric vector 
##' @param X2 numeric vector 
##' @param l length parameter
##' 
##' @return the covariance matrix of X1 and X2
Exp_k <- function(X1,X2,l=1){
  Sigma <- matrix(rep(0, length(X1)*length(X2)), nrow=length(X1))
  for (i in 1:nrow(Sigma)) {
    for (j in 1:ncol(Sigma)) {
      r <- abs(X1[i]-X2[j])
      Sigma[i,j] <- exp(-r/l)
    }
  }
  return (Sigma)
}

##' 3/2 Matern kernel 
##' 
##' @param X1 numeric vector 
##' @param X2 numeric vector 
##' @param l length parameter
##' 
##' @return the covariance matrix of X1 and X2
Matern_k <- function(X1,X2,l=1){
  Sigma <- matrix(rep(0, length(X1)*length(X2)), nrow=length(X1))
  for (i in 1:nrow(Sigma)) {
    for (j in 1:ncol(Sigma)) {
      r <- abs(X1[i]-X2[j])
      Sigma[i,j] <- (1+sqrt(3)*r/l)*exp(-sqrt(3)*r/l)
    }
  }
  return (Sigma)
}

##' The objective function f constructed by 3/2 Matern
##' 
##' @param x
##' 
##' @return f(x)
f <- function(x,r){
  l <- length(x)
  Sigma <- Matern_k(x,x,r)
  R <- chol(Sigma)
  return( t(R)%*%c(rnorm(n=l,mean=0,sd = 1)) )
}


##' The objective function g constructed by Exp
##' 
##' @param x
##' 
##' @return g(x)
g <- function(x,r){
  l <- length(x)
  Sigma <- Exp_k(x,x,r)
  R <- chol(Sigma)
  return( t(R)%*%c(rnorm(n=l,mean=0,sd = 1))  )
}

