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

##' Calculate the predictive probability (!Does Not Work!)
##' probability measure for the partitions
##' 
##' @param alpha vector of parameter alpha
##' @param maxK scalar, maximum tree levels i.e.depth of the tree
##' 
##' @return vector of expected predictive probabilities 
##' for each partition
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
    alpha[2*m-1] <- m^2
    alpha[2*m] <- m^2
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
  log_marg_prev <- c(-Inf)
  log_marg_after <- c(-Inf)
  for (tau in 2:n_obs){
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
