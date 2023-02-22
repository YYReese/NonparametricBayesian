source("functions.R")

# observations
n_obs <- 500
t <- 1:n_obs

# change in mean in Exp
y_exp1 <- rexp(n_obs/2) 
y_exp2 <- rexp(n_obs/2)+1
obs_exp <- c(y_exp1,y_exp2)
# change in mean in Gaussian
y_gaussian1 <- rnorm(n_obs/2,mean=0,sd=.5)
y_gaussian2 <- rnorm(n_obs/2,mean=2,sd=.5)
obs_gaussian_mean <- c(y_gaussian1,y_gaussian2)
# change in variance in Gaussian
y_gaussian3 <- rnorm(n_obs/2,mean=0,sd=2)
obs_gaussian_var <- c(y_gaussian1,y_gaussian3)

# initialisation
maxK <- 7 # maximum depth
a_prior <- alpha_priors(maxK=maxK) # alpha prior


##############Single change in mean --Exp##########################
tau_exp <- find_change_location(obs_exp, a_prior=a_prior)
tau_gaussian_mean <- find_change_location(obs_gaussian_mean,
                                          a_prior=a_prior)
tau_gaussian_var <- find_change_location(obs_gaussian_var,
                                         a_prior=a_prior)

par(mfrow=c(1,3))
plot(obs_exp,ylab="Observed value",
     main="Change in mean in exponential data")
abline(v=tau_exp,col="blue")
abline(v=n_obs/2,col="red")

plot(obs_gaussian_mean,ylab="Observed value",
     main="Change in mean in Gaussian data")
abline(v=tau_gaussian_mean,col="blue")
abline(v=n_obs/2,col="red")

plot(obs_gaussian_var,ylab="Observed value",
     main="Change in variance in Gaussian data")
abline(v=tau_gaussian_var,col="blue")
abline(v=n_obs/2,col="red")


#############M1:Single change point detection--non-Gaussian############
##################Work well as desired##############################
library(readr)
machine_0 <- read_csv("data/machine_0.csv")

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
points(change1,log_marg[change1],col="red")


##################Another example#################
library("rjson")
brent_spot <- fromJSON(file="brent_spot.json")
brent_spot <- as.data.frame(brent_spot)
plot(brent_spot$series.raw)

# observations
n_obs <- nrow(stationB)
t <- 1:n_obs
y_obs <- stationB$PRCP

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
points(change1,log_marg[change1],col="red")

plot(stationB$PRCP,type="l")
abline(v=change1, col="red")


