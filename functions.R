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

##' Calculate the updated/posterior probability measure for the partitions
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

##' Draw samples from the posterior of theta
##' 
##' @param a_post vector of alpha_post
##' 
##' @return vector of theta drawing from the posterior 
draw_theta_post <- function(a_post){
  theta <- c()
  for(i in 1:(length(a_post)/2)){
    theta <- c(theta, rbeta(1, a_post[2*i-1], a_post[2*i]))
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
  theta <- draw_theta_post(a_post)
  
  prob <- Cal_prob(theta, maxK)
  
}
