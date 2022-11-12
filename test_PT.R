# observations
y_obs <- rnorm(100000)
# Interval
a <- -5
b <- 5

# maximum depth of the tree
maxLevel <- 10

# Break the whole set into disjoint partitions
#(each as a tree node)
breaks<-a+(1:(2^maxLevel-1))*(b-a)/(2^maxLevel)

# Count the number of observations falling into each partition
n_counts <- hist(y_obs,breaks=c(-Inf,breaks,Inf),plot=FALSE)$counts
number_of_nodes <- length(n_counts)

# Set the prior knowledge about alpha_l and alpha_r
# Default: vector of 1s
a_l_prior <- rep(1/2,number_of_nodes)
a_r_prior <- rep(1/2, number_of_nodes)

# Update the posterior of alpha_l and alpha_r
a_l <- a_l_prior
a_r <- a_r_prior
for (yi in y_obs){
  a_l <- update_alpha(a_l,a_r,maxLevel,yi)$a_l_post
  a_r <- update_alpha(a_l,a_r,maxLevel,yi)$a_r_post
}

# Draw samples from the posterior of theta
theta <- rep(NA,number_of_nodes)
for (i in 1:number_of_nodes){
  theta[i] <- rbeta(1,a_l[i],a_r[i])
}


# Calculate the probability measure for a new y
Cal_F(a,b,theta, maxLevel, breaks, 0)

# n=1000 new y
y_new <- seq(-3,3, len=100)
p <- rep(1,length(y_new))
for (i in 1:length(y_new)){
  p[i] <- Cal_F(a,b,theta,maxLevel, breaks,y_new[i])
}


ggplot(data.frame(x=y_new,y=p/sum(p),y_true=dnorm(y_new))) +
  geom_line(aes(x,y)) +
  geom_line(aes(x,y_true),col="red")


################################################################################

# Calculate the marginal likelihood of y_observation
# Use eqn(6) 
marginal_likeli(y=y, a_l=a_l_prior,a_r=a_r_prior,n_counts = n_counts)


prob <- rep(NA, length(breaks))
for (i in length(breaks)){
  prob[i] <- prob_measure(a_l_post,a_r_post,breaks,breaks[i])
}

plot(breaks,prob)
# Calculate the posterior density of 
post_density <- function(p_measure,likeli){
  p_measure/likeli
}
