source("functions.R")

x <- 1:3000
dx <- (x[length(x)]-x[1])/length(x)
fx <- f(x)
#gx <- g(x)
y1 <- fx+rnorm(length(x),sd=0.1)

#y2 <- gx+rnorm(length(x),sd=0.1)

y <- y1
# maximum depth of the tree
maxK <- 7

# Interval
a <- floor(min(y))
b <- ceiling(max(y))

nParam <- 2^(maxK + 1) - 2
breaks <- a + (1:(2^maxK-1)) / (2^maxK) * (b-a) 
# Set the prior knowledge about alpha
a_prior <- rep(1/2,  nParam)

# observations
n_obs <- 2000
obs_idx <- runif(n_obs,min=0,max=length(y))
y_obs <- y[obs_idx]


# Draw n realisations from posterior PT
n <- 1000
prob <- matrix(nrow=n,ncol=nParam)
for (i in 1:n){
  prob[i,] <- PT_update(a,b,maxK,a_prior, y_obs)$prob
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
     xlab="y",xlim=c(min(y),max(y)),
     ylim=c(0, max(up_p/((b-a)/2^maxK))), ylab="Density")
# 90% Confidence interval
lines(xPoints, up_p/((b-a)/2^maxK),col="grey",type="l")
lines(xPoints, lw_p/((b-a)/2^maxK),col="grey",type="l")


# Show true density
lines(density(y), col=3,pch=3, cex=0.05)


title("n=5000 observations",cex=0.3)

legend("topright", pch=c(15,15), 
       legend=c("True density","Observations" ,"Post PT estimation","90% CI"), 
       col=c(3,4,2,"grey"), bty="n", cex=0.5)

