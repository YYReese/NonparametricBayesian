# Interval
a <- 0
b <- 8

# maximum depth of the tree
maxK <- 10

nParam <- 2^(maxK + 1) - 2
breaks <- a + (1:(2^maxK-1)) / (2^maxK) * (b-a) 

# Set the prior knowledge about alpha
a_prior <- rep(1/2,  nParam)

# observations
y_obs <- rpois(10000, 2)

# Draw n realisations from PT 
n <- 1000
prob <- matrix(nrow=n,ncol=nParam)
for (i in 1:n){
  prob[i,] <- PT_update(a,b,maxK,a_prior, y_obs)
}

# Compute the mean of the realisations at each point
mu_p <- rep(NA, ncol(prob))
for (j in 1:ncol(prob)){
  mu_p[j] <- mean(prob[,j])
}

xPoints <- (c(a, breaks) + c(breaks, b)) / 2
p <- mu_p[(length(mu_p)-2^(maxK)+1): length(mu_p)]

# Plot mean of the n realisations from PT
plot(xPoints, p, col=2, type="p", 
     ylim=c(-0.7, 0.7), ylab="Density")

# Plot some of the realisations
#for (i in 1:2){
#  points(xPoints, prob[i,(length(mu_p)-2^(maxK)+1): length(mu_p)]/((b-a)/2^maxK),type="p", 
#        ylim=c(-0.7, 0.7), ylab="Density")
#}

# Show the true distribution 
points(a:b, dpois(a:b,2), col=3, type="l")
legend("topright", pch=c(15,15), 
       legend=c("PT Tree","Poisson"), col=c(2,3), bty="n")
