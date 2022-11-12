##' Calculate the probability measure given theta
##' 
##' @param y new data point
##' @param theta vector of parameter theta
##' @param maxLevel scalar, maximum tree levels i.e.depth of the tree
##' @param breaks vector of break points
##' 
##' @return scalar, probability measure of y
Cal_F <- function(a, b, theta, maxLevel, breaks, y){
  Mid <- (length(breaks)+1)/2
  p <- theta[Mid]
  for (i in 2:maxLevel-1) {
    delta <- (b-a)/2^(i+1)
    if (y < breaks[Mid]){
      Mid <- Mid - delta
      p <- p * theta[Mid]
    }
    else {
      Mid <- Mid + delta
      p <- p * (1-theta[Mid])
    }
  }
  return (p)
}

##' Find path
##' @param m depth of the tree
##' @param breaks vector of break points
##' @param y new data point
Find_path <- function(a,b,m, breaks, y){
  idx <- rep(FALSE, length(breaks))
  Mid <- (length(breaks)+1)/2
  idx[Mid] <- TRUE
  for (i in 1:m) {
    delta <- (b-a)/2^(i+1)
    if (y < breaks[Mid]){
      Mid <- Mid - delta
      idx[Mid] <- TRUE
    }
    else {
      Mid <- Mid + delta
      idx[Mid] <- TRUE
    }
  }
  return (idx)
}


##' Update alpha
##' @param a_l vector of scalar, alpha_l prior
##' @param a_r vector of scalar, alpha_r prior
##' @param y single observation
##' @param m the depth of the tree
##' 
##' @return a_l_post, a_r_post
update_alpha <- function(a_l,a_r,m,y){
  Mid <- (length(breaks)+1)/2
  for (i in 1:(m-1)) {
    delta <- 2^(m-i-1)
    if (y < breaks[Mid]){
      a_l[Mid] <- a_l[Mid]+1
      Mid <- Mid - delta
    }
    else {
      a_r[Mid] <- a_r[Mid]+1
      Mid <- Mid + delta
    }
  }
  return (list(a_l_post=a_l,a_r_post=a_r))
}

update_alpha(a_l_prior,a_r_prior,maxLevel,0.13)

































##' Calculate the marginal likelihood for y
##' 
##' @param y vector of observations
##' @param a_l vector of prior values of alpha_epsilon_k_0 
##' (parameter 1 for theta in Beta distribution)
##' note: theta ~ Beta(a_l,a_r)
##' @param a_r vector of prior values of alpha_epsilon_k_1
##' (parameter 2 for theta in Beta distribution)
##' @param n_counts vector of number of observations in each node/sub_interval
##' 
##' @return scalar, the marginal likelihood of observations y
marginal_likeli <- function(y, a_l, a_r, n_counts){
  log_likeli <- lgamma(a_l+a_r) + lgamma(a_l+n_counts) + lgamma(a_r+n_counts) - 
    lgamma(a_l) - lgamma(a_r) - lgamma(a_r+a_r+2*n_counts)
  return (exp(log_likeli))
}

##' Calculate the probability measure of each Bk
##' 
##' @param a_l 
##' @param a_r
##' @param breaks "tree nodes"
##' @param y scalar, new data point where we would like to measure the probability
##' 
##' @return probability measure of y 
prob_measure <- function(a_l, a_r, breaks, y){
  prod(a_l[y<tail(breaks,length(breaks)-1)]/
         (a_l[y<tail(breaks,length(breaks)-1)]+
            a_r[y<tail(breaks,length(breaks)-1)])) * 
    prod(a_r[y>=tail(breaks,length(breaks)-1)]/
           (a_l[y>=tail(breaks,length(breaks)-1)]+
              a_r[y>=tail(breaks,length(breaks)-1)]))
}