# Interval
a <- -5
b <- 5



# maximum depth of the tree
maxK <- 10
nParam <- 2^(maxK + 1) - 2

# Break the whole set into disjoint partitions (leaf)
breaks <- a + (1:(2^maxK-1)) / (2^maxK) * (b-a) 


# Set the prior knowledge about alpha
a_prior <- rep(1/2, nParam)




# observations
y_obs <- rnorm(1000000)
y_obs[y_obs < a] = a
y_obs[y_obs > b] = b

# Count the number of observations falling into each partition
leaf_count <- hist(y_obs, breaks=c(-Inf, breaks, Inf), plot=FALSE)$counts


# find the carent count
get_parent_count <- function(child_count){
  parent_len = length(child_count) / 2
  return (child_count[2*(1:parent_len)-1] + child_count[2*(1:parent_len)])
}


count = c()
child_count = leaf_count
while(length(child_count) > 1){
  count = c(child_count, count)
  parent_count = get_parent_count(child_count)
  child_count = parent_count
}
   
  
a_post = a_prior + count



# Draw samples from the posterior of theta
theta = c()
for(i in 1:(length(a_post)/2)){
  theta = c(theta, rbeta(1, a_post[2*i-1], a_post[2*i]))
}


prob = c(theta[1], 1-theta[1])
prob_i = c(theta[1], 1-theta[1])
for(i in 2:maxK){
  new_prob_i = c()
  for(j in 1:length(prob_i)){
    x = theta[2^(i-1) + j - 1]
    new_prob_i = c(new_prob_i, prob_i[j]*x, prob_i[j]*(1-x))
  }
  prob_i = new_prob_i
  prob = c(prob, new_prob_i)
}




# compare density
get_prob <- function(x, breaks, leaf_prob){
  i = 1
  iEnd = length(leaf_prob)
  while((x > breaks[i]) & (i < iEnd)){
    i = i + 1
  }
  return(leaf_prob[i])
}


xPoints <- (c(a, breaks) + c(breaks, b)) / 2
p = prob[(length(prob)-2^(maxK)+1): length(prob)]

plot(xPoints, p/((b-a)/2^maxK), col=2, type="l", ylim=c(-0.7, 0.7), ylab="Density")
lines(xPoints, dnorm(xPoints), col=3, type="l")
legend("topright", pch=c(15,15), legend=c("PT Tree","Norm"), col=c(2,3), bty="n")

