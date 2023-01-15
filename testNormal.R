source("functions.R")
# Interval
a <- -5
b <- 5

# maximum depth of the tree
maxK <- 7

nParam <- 2^(maxK + 1) - 2
breaks <- a + (1:(2^maxK-1)) / (2^maxK) * (b-a) 
# Set the prior knowledge about alpha
a_prior <- rep(1/2,  nParam)


# observations
n_obs <- 1000
y_obs <- rnorm(n_obs) 

# Draw n realisations from posterior PT and prior PT
n <- 10000
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
plot(xPoints, p/((b-a)/2^maxK), col=2, type="l", 
     xlab="y",
     ylim=c(-0.7, 0.7), ylab="Density")

# 90% Confidence interval
lines(xPoints, up_p/((b-a)/2^maxK),col="grey",type="l")
lines(xPoints, lw_p/((b-a)/2^maxK),col="grey",type="l")

# Show the observation points 
points(density(y_obs), col=4,pch=3, cex=0.05)

# Standard Gaussian
lines(xPoints,dnorm(xPoints))


title("n=1000 observations",cex=0.3)

legend("topright", pch=c(15,15), 
       legend=c("Prior PT", "Post PT","Observations","90% CI"), 
       col=c("black",2,4,"grey"), bty="n", cex=0.5)

