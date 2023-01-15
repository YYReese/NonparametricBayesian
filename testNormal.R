source("functions.R")
# Interval
a <- -5
b <- 5

# maximum depth of the tree
maxK <- 5

nParam <- 2^(maxK + 1) - 2
breaks <- a + (1:(2^maxK-1)) / (2^maxK) * (b-a) 
# Set the prior knowledge about alpha
a_prior <- rep(1/2,  nParam)

# observations
n_obs <- 1000
y_obs <- rnorm(n_obs)

# Draw n realisations from posterior PT 
n <- 1000
prob <- matrix(nrow=n,ncol=nParam)
for (i in 1:n){
  prob[i,] <- PT_update(a,b,maxK,a_prior, y_obs)
}

# Compute the mean of the realisations at each point
mu_p <- rep(NA, ncol(prob))
up_p <- rep(NA, ncol(prob))
lw_p <- rep(NA, ncol(prob))
for (j in 1:ncol(prob)){
  mu_p[j] <- mean(prob[,j])
  up_p[j] <- quantile(prob[,j],prob=0.95,type=1)
  lw_p[j] <- quantile(prob[,j],prob=0.05,type=1)
}

xPoints <- (c(a, breaks) + c(breaks, b)) / 2
p <- mu_p[(length(mu_p)-2^(maxK)+1): length(mu_p)]
up_p <- up_p[(length(up_p)-2^(maxK)+1): length(up_p)]
lw_p <- lw_p[(length(lw_p)-2^(maxK)+1): length(lw_p)]

# Plot mean of the n realizations from PT
plot(xPoints, p/((b-a)/2^maxK), col=2, type="s", 
     ylim=c(-0.7, 0.7), ylab="Density")

# Plot a sample from prior PT
lines(xPoints, Cal_prob(draw_theta(a_prior), maxK)
     [(length(mu_p)-2^(maxK)+1): length(mu_p)],type = "s")

# 90% Confidence interval
lines(xPoints, up_p/((b-a)/2^maxK),col="grey")
lines(xPoints, lw_p/((b-a)/2^maxK),col="grey")

# Show the true data distribution 
lines(density(y_obs), col=3, type="l")

legend("topright", pch=c(15,15), 
       legend=c("PT","True density","90% CI"), 
       col=c(2,3,"grey"), bty="n", cex=0.5)
