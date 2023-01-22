source("functions.R")
load("simulatedData.RData")

# observations
n_obs <- 5000
t <- 1:n_obs
y_obs1 <- rnorm(n_obs/2) 
y_obs2 <- rexp(n_obs/2)
y_obs <- c(y_obs1,y_obs2)

# maximum depth of the tree
maxK <- 7

# Interval
a <- floor(min(y_obs))
b <- ceiling(max(y_obs))

nParam <- 2^(maxK + 1) - 2
breaks <- a + (1:(2^maxK-1)) / (2^maxK) * (b-a) 
# Set the prior knowledge about alpha
a_prior <- rep(1/2,  nParam)

#########Method 1:Single change point detection--Gaussian##################
###############Work well as desired but ############################
log_marg_prev <- c(-Inf)
log_marg_after <- c(-Inf)
for (tau in 2:n_obs){
  y_obs_prev <- y_obs[1:tau-1]
  y_obs_after <- y_obs[tau:n_obs]

  a_post_prev <- PT_update(a,b,maxK,a_prior, y_obs_prev)$a_post

  a_post_after <- PT_update(a,b,maxK,a_prior, y_obs_after)$a_post
  log_marg_prev[tau] <- log_Marginal_prob(a_prior,a_post_prev)
  log_marg_after[tau] <- log_Marginal_prob(a_prior,a_post_after)
}
log_marg <- log_marg_prev + log_marg_after
change <- which.max(log_marg)
plot(log_marg)
points(change,log_marg[change],col="red")

#########Method 2:Single change point detection--Gaussian##################
##################Does not work but makes more sense################
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
plot(log_marg)



#############M1:Single change point detection--non-Gaussian############
##################Work well as desired##############################
library(readr)
machine_0 <- read_csv("machine_0.csv")

# observations
n_obs <- 3000
t <- 1:n_obs
y_obs <- machine_0$`3`

# maximum depth of the tree
maxK <- 7

# Interval
a <- floor(min(y_obs))
b <- ceiling(max(y_obs))

nParam <- 2^(maxK + 1) - 2
breaks <- a + (1:(2^maxK-1)) / (2^maxK) * (b-a) 
# Set the prior knowledge about alpha
a_prior <- rep(1/2,  nParam)
log_marg_prev <- c(-Inf)
log_marg_after <- c(-Inf)
for (tau in 2:n_obs){
  # Observations
  y_obs_prev <- y_obs[1:tau-1]
  y_obs_after <- y_obs[tau:n_obs]
  
  # Interval
  a <- floor(min(y_obs_prev))
  b <- ceiling(max(y_obs_prev))
  # Breaks
  nParam <- 2^(maxK + 1) - 2
  breaks <- a + (1:(2^maxK-1)) / (2^maxK) * (b-a) 
  # Fitting
  a_post_prev <- PT_update(a,b,maxK,a_prior, y_obs_prev)$a_post
  # Interval
  a <- floor(min(y_obs_after))
  b <- ceiling(max(y_obs_after))
  # Breaks
  breaks <- a + (1:(2^maxK-1)) / (2^maxK) * (b-a) 
  a_post_after <- PT_update(a,b,maxK,a_prior, y_obs_after)$a_post
  log_marg_prev[tau] <- log_Marginal_prob(a_prior,a_post_prev)
  log_marg_after[tau] <- log_Marginal_prob(a_prior,a_post_after)
}
log_marg <- log_marg_prev + log_marg_after
change1 <- which.max(log_marg)
plot(log_marg)
points(change1,log_marg[change],col="red")



#############M2:Single change point detection--non-Gaussian############
##################Work well as desired##############################
# observations
n_obs <- 3000
t <- 1:n_obs
y_obs <- machine_0$`3`

# maximum depth of the tree
maxK <- 7

# Interval
a <- floor(min(y_obs))
b <- ceiling(max(y_obs))
a_prior <- rep(1/2,  nParam)
log_marg_prev <- c(-Inf)
log_marg_after <- c(-Inf)
for (tau in 2:n_obs){
  y_obs_prev <- y_obs[1:tau-1]
  y_obs_after <- y_obs[tau:n_obs]
  
  a_post_prev <- PT_update(a,b,maxK,a_prior, y_obs_prev)$a_post
  
  a_post_after <- PT_update(a,b,maxK,a_prior, y_obs_after)$a_post
  log_marg_prev[tau] <- log_Marginal_prob(a_prior,a_post_prev)
  log_marg_after[tau] <- log_Marginal_prob(a_prior,a_post_after)
}
log_marg <- log_marg_prev + log_marg_after
change2 <- which.max(log_marg)
plot(log_marg)
points(change2,log_marg[change],col="red")

plot(machine_0$`3`)
abline(v=change1, col="red")
abline(v=change1, col="blue")




