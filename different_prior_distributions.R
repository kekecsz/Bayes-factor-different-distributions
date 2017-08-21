
##################################################
#                                                #
#                Custom functions                #
#                                                #
##################################################

## these are functions calculating Bayes factor using different prior distributions
## in all of these functions information can be provided about results of a previous study
## to be able to calculate replication Bayes factor. Or uninformed prior can also be used.


### Functions for beta prior

fullAlt_beta = Vectorize(function(p, y, N, alpha, beta){
  exp(dbinom(y, N, p, log = TRUE) + dbeta(p, alpha, beta, log = TRUE)) 
},"p")

normalize_beta = function(alpha, beta, interval){
  diff(pbeta(interval, alpha, beta))
}

restrictedAlt_beta = function(p,y,N,y_prior,N_prior,interval){
  alpha = y_prior + 1
  beta = N_prior - y_prior + 1
  fullAlt_beta(p, y, N, alpha, beta) / normalize_beta(alpha, beta, interval) * (p>interval[1] & p<interval[2])
}

margLike_beta = function(y, N, y_prior, N_prior, interval){
  integrate(restrictedAlt_beta, interval[1], interval[2], 
            y = y, N = N, y_prior = y_prior, N_prior = N_prior, interval = interval)[[1]]
}

BF01_beta = Vectorize(function(y, N, y_prior, N_prior, interval, null_prob){
  dbinom(y, N, null_prob) / margLike_beta(y = y, N = N, y_prior = y_prior, N_prior = N_prior, interval = interval)
},"y")






### Funtion for normal prior approximating the beta distribution

fullAlt_norm_p = Vectorize(function(p, y, N, mean, sd){
  exp(dbinom(y, N, p, log = TRUE) + dnorm(p, mean, sd, log = TRUE)) 
},"p")

normalize_norm_p = function(mean, sd, interval){
  diff(pnorm(interval, mean, sd))
}

restrictedAlt_norm_p = function(p,y,N,y_prior,N_prior,interval){
  alpha = y_prior + 1
  beta = N_prior - y_prior + 1
  mean = alpha/(alpha + beta)
  var = (alpha*beta)/(((alpha + beta)^2)*(1 + alpha + beta))
  sd = sqrt(var)
  
  fullAlt_norm_p(p, y, N, mean, sd) / normalize_norm_p(mean, sd, interval) * (p>interval[1] & p<interval[2])
}

margLike_norm_p = function(y, N, y_prior, N_prior, interval){
  integrate(restrictedAlt_norm_p, interval[1], interval[2], 
            y = y, N = N, y_prior = y_prior, N_prior = N_prior, interval = interval)[[1]]
}

BF01_norm_approx_beta = Vectorize(function(y, N, y_prior, N_prior, interval, null_prob){
  dbinom(y, N, null_prob) / margLike_norm_p(y = y, N = N, y_prior = y_prior, N_prior = N_prior, interval = interval)
},"y")



### functions for logistic prior
### based on Richard Morey's functions

fullAlt_logis = Vectorize(function(theta, y, N, prior_mean_lo, prior_sd_lo){
  p = plogis(theta)
  exp(dbinom(y, N, p, log = TRUE) + dlogis(theta, prior_mean_lo, prior_sd_lo, log = TRUE)) # here, scale is expected in cohen's d which is converted to log odds. theta and shift is already in log odds, shift is the mean of the normal prior (in log odds). 
},"theta")

normalize_logis = function(prior_mean_lo, prior_sd_lo, interval){
  diff(plogis(interval, prior_mean_lo, prior_sd_lo))
}

restrictedAlt_logis = function(theta,y,N,y_prior,N_prior,interval){
  alpha = y_prior + 1
  beta = N_prior - y_prior + 1
  mean = alpha/(alpha + beta)
  var = (alpha*beta)/(((alpha + beta)^2)*(1 + alpha + beta))
  prior_mean_lo = qlogis(mean)
  prior_sd_lo = sqrt( var / dlogis( prior_mean_lo )^2 )
  
  fullAlt_logis(theta,y,N,prior_mean_lo,prior_sd_lo) / normalize_logis(prior_mean_lo, prior_sd_lo, interval) * (theta>interval[1] & theta<interval[2])
}

margLike_logis = function(y, N, y_prior, N_prior, interval){
  theta_interval = qlogis(sort(interval))
  integrate(restrictedAlt_logis, theta_interval[1], theta_interval[2], 
            y = y, N = N, y_prior = y_prior, N_prior = N_prior, interval = theta_interval)[[1]]
}

BF01_logis = Vectorize(function(y, N, y_prior, N_prior, interval, null_prob){
  dbinom(y,N,null_prob) / margLike_logis(y, N, y_prior, N_prior, interval)
},"y")




### functions for normal prior approximating the logistic distibution
### based on Richard Morey's functions

fullAlt_norm = Vectorize(function(theta, y, N, prior_mean_lo, prior_sd_lo){
  p = plogis(theta)
  exp(dbinom(y, N, p, log = TRUE) + dnorm(theta, prior_mean_lo, prior_sd_lo, log = TRUE)) # here, scale is expected in cohen's d which is converted to log odds. theta and shift is already in log odds, shift is the mean of the normal prior (in log odds). 
},"theta")

normalize_norm = function(prior_mean_lo, prior_sd_lo, interval){
  diff(pnorm(interval, prior_mean_lo, prior_sd_lo))
}

restrictedAlt_norm = function(theta,y,N,y_prior,N_prior,interval){
  alpha = y_prior + 1
  beta = N_prior - y_prior + 1
  mean = alpha/(alpha + beta)
  var = (alpha*beta)/(((alpha + beta)^2)*(1 + alpha + beta))
  prior_mean_lo = qlogis(mean)
  prior_sd_lo = sqrt( var / dlogis( prior_mean_lo )^2 )* pi/sqrt(3)
  
  fullAlt_norm(theta,y,N,prior_mean_lo,prior_sd_lo) / normalize_norm(prior_mean_lo, prior_sd_lo, interval) * (theta>interval[1] & theta<interval[2])
}

margLike_norm = function(y, N, y_prior, N_prior, interval){
  theta_interval = qlogis(sort(interval))
  integrate(restrictedAlt_norm, theta_interval[1], theta_interval[2], 
            y = y, N = N, y_prior = y_prior, N_prior = N_prior, interval = theta_interval)[[1]]
}

BF01_norm_approx_logis = Vectorize(function(y, N, y_prior, N_prior, interval, null_prob){
  dbinom(y,N,null_prob) / margLike_norm(y, N, y_prior, N_prior, interval)
},"y")








##################################################
#                                                #
#                 Demonstrations                 #
#                                                #
##################################################



##### demonstration with small positive effect in the new study and uninformed prior with the same center as H0

observed_p = 0.51

multip = 1:300 # to simulate different sample sizes with the same observed p
ys <- 100*observed_p*multip
Ns <- 100*multip

bf.beta = NA
bf.logis = NA
bf.norm_approx_logis = NA
bf.norm_approx_beta = NA

null_prob = 0.5
y_prior = 6.5 # information from a previous study. To get uninformed prior: y_prior = 6.5 and N_prior = 13 will result in rscale = 0.5, y_prior = 0.5 and N_prior = 1 will result in rscale = 1 in the proportionBF function of the BayesFactor package
N_prior = 13 # information from a previous study. To get uninformed prior: y_prior = 6.5 and N_prior = 13 will result in rscale = 0.5, y_prior = 0.5 and N_prior = 1 will result in rscale = 1 in the proportionBF function of the BayesFactor package
interval = c(0.5,1)

# the below is needed to calculate the rscale parameter for the proportionBF within BayesFactor package
alpha = y_prior + 1
beta = N_prior - y_prior + 1
mean = alpha/(alpha + beta)
var = (alpha*beta)/(((alpha + beta)^2)*(1 + alpha + beta))
prior_mean_lo = qlogis(mean)
prior_sd_lo = sqrt( var / dlogis( prior_mean_lo )^2 )

for(i in 1:length(multip)){
  y = ys[i]
  N = Ns[i]
  
  bf.beta[i] = 1/BF01_beta(y = y, N = N, y_prior = y_prior, N_prior = N_prior, interval = interval, null_prob = null_prob)
  bf.logis[i] = 1/BF01_logis(y = y, N = N, y_prior = y_prior, N_prior = N_prior, interval = interval, null_prob = null_prob)
  bf.norm_approx_logis[i] = 1/BF01_norm_approx_logis(y = y, N = N, y_prior = y_prior, N_prior = N_prior, interval = interval, null_prob = null_prob)
  bf.norm_approx_beta[i] = 1/BF01_norm_approx_beta(y = y, N = N, y_prior = y_prior, N_prior = N_prior, interval = interval, null_prob = null_prob)
 
}  

par( mfrow = c( 1, 2 ) )

# plot for small sample sizes
plot(ys[1:20],bf.beta[1:20],ty='n',log="y",ylab="Bayes factor",xlab="y",las=1, ylim = c(min(c(bf.beta[1:20], bf.logis[1:20], bf.norm_approx_logis[1:20], bf.norm_approx_beta[1:20])),max(c(  bf.beta[1:20], bf.logis[1:20], bf.norm_approx_logis[1:20], bf.norm_approx_beta[1:20]))))
lines(ys[1:20],bf.beta[1:20],col="black") # BF from the beta distribution
lines(ys[1:20],bf.logis[1:20],col="blue") # normal approcimation of the beta distribution, should match BF from beta distribution if y and N and prior y and N are high enough
lines(ys[1:20],bf.norm_approx_beta[1:20],col="green") # BF from the logistic distribution, same as used in the BayesFactor package in the proportionBF function
lines(ys[1:20],bf.norm_approx_logis[1:20],col="red") # normal approximation of the logistic prior used in the BayesFactor package in the proportionBF function, should match the BFs from logistic distributions closely, but it is not a perfect match

# plot for large sample sizes
plot(ys[290:300],bf.beta[290:300],ty='n',log="y",ylab="Bayes factor",xlab="y",las=1, ylim = c(min(c(bf.beta[290:300], bf.logis[290:300], bf.norm_approx_logis[290:300], bf.norm_approx_beta[290:300])),max(c(  bf.beta[290:300], bf.logis[290:300], bf.norm_approx_logis[290:300], bf.norm_approx_beta[290:300]))))
lines(ys[290:300],bf.beta[290:300],col="black") # BF from the beta distribution
lines(ys[290:300],bf.logis[290:300],col="blue") # normal approcimation of the beta distribution, should match BF from beta distribution if y and N and prior y and N are high enough
lines(ys[290:300],bf.norm_approx_beta[290:300],col="green") # BF from the logistic distribution, same as used in the BayesFactor package in the proportionBF function
lines(ys[290:300],bf.norm_approx_logis[290:300],col="red") # normal approximation of the logistic prior used in the BayesFactor package in the proportionBF function, should match the BFs from logistic distributions closely, but it is not a perfect match







##### demonstration with null effect in the new study and uninformed prior with the same center as H0

observed_p = 0.50

multip = 1:300 # to simulate different sample sizes with the same observed p
ys <- 100*observed_p*multip
Ns <- 100*multip

bf.beta = NA
bf.logis = NA
bf.norm_approx_logis = NA
bf.norm_approx_beta = NA

null_prob = 0.5
y_prior = 6.5 # information from a previous study. To get uninformed prior: y_prior = 6.5 and N_prior = 13 will result in rscale = 0.5, y_prior = 0.5 and N_prior = 1 will result in rscale = 1 in the proportionBF function of the BayesFactor package
N_prior = 13 # information from a previous study. To get uninformed prior: y_prior = 6.5 and N_prior = 13 will result in rscale = 0.5, y_prior = 0.5 and N_prior = 1 will result in rscale = 1 in the proportionBF function of the BayesFactor package
interval = c(0.5,1)

# the below is needed to calculate the rscale parameter for the proportionBF within BayesFactor package
alpha = y_prior + 1
beta = N_prior - y_prior + 1
mean = alpha/(alpha + beta)
var = (alpha*beta)/(((alpha + beta)^2)*(1 + alpha + beta))
prior_mean_lo = qlogis(mean)
prior_sd_lo = sqrt( var / dlogis( prior_mean_lo )^2 )

for(i in 1:length(multip)){
  y = ys[i]
  N = Ns[i]
  
  bf.beta[i] = 1/BF01_beta(y = y, N = N, y_prior = y_prior, N_prior = N_prior, interval = interval, null_prob = null_prob)
  bf.logis[i] = 1/BF01_logis(y = y, N = N, y_prior = y_prior, N_prior = N_prior, interval = interval, null_prob = null_prob)
  bf.norm_approx_logis[i] = 1/BF01_norm_approx_logis(y = y, N = N, y_prior = y_prior, N_prior = N_prior, interval = interval, null_prob = null_prob)
  bf.norm_approx_beta[i] = 1/BF01_norm_approx_beta(y = y, N = N, y_prior = y_prior, N_prior = N_prior, interval = interval, null_prob = null_prob)

}  

par( mfrow = c( 1, 2 ) )

# plot for small sample sizes
plot(ys[1:20],bf.beta[1:20],ty='n',log="y",ylab="Bayes factor",xlab="y",las=1, ylim = c(min(c(bf.beta[1:20], bf.logis[1:20], bf.norm_approx_logis[1:20], bf.norm_approx_beta[1:20])),max(c(  bf.beta[1:20], bf.logis[1:20], bf.norm_approx_logis[1:20], bf.norm_approx_beta[1:20]))))
lines(ys[1:20],bf.beta[1:20],col="black") # BF from the beta distribution
lines(ys[1:20],bf.logis[1:20],col="blue") # normal approcimation of the beta distribution, should match BF from beta distribution if y and N and prior y and N are high enough
lines(ys[1:20],bf.norm_approx_beta[1:20],col="green") # BF from the logistic distribution, same as used in the BayesFactor package in the proportionBF function
lines(ys[1:20],bf.norm_approx_logis[1:20],col="red") # normal approximation of the logistic prior used in the BayesFactor package in the proportionBF function, should match the BFs from logistic distributions closely, but it is not a perfect match

# plot for large sample sizes
plot(ys[290:300],bf.beta[290:300],ty='n',log="y",ylab="Bayes factor",xlab="y",las=1, ylim = c(min(c(bf.beta[290:300], bf.logis[290:300], bf.norm_approx_logis[290:300], bf.norm_approx_beta[290:300])),max(c(  bf.beta[290:300], bf.logis[290:300], bf.norm_approx_logis[290:300], bf.norm_approx_beta[290:300]))))
lines(ys[290:300],bf.beta[290:300],col="black") # BF from the beta distribution
lines(ys[290:300],bf.logis[290:300],col="blue") # normal approcimation of the beta distribution, should match BF from beta distribution if y and N and prior y and N are high enough
lines(ys[290:300],bf.norm_approx_beta[290:300],col="green") # BF from the logistic distribution, same as used in the BayesFactor package in the proportionBF function
lines(ys[290:300],bf.norm_approx_logis[290:300],col="red") # normal approximation of the logistic prior used in the BayesFactor package in the proportionBF function, should match the BFs from logistic distributions closely, but it is not a perfect match








##### demonstration with small positive effect in the new study and prior distribution based on results of a previous study finding a small positive effect
observed_p = 0.51

multip = 1:300 # to simulate different sample sizes with the same observed p
ys <- 100*observed_p*multip
Ns <- 100*multip

bf.beta = NA
bf.logis = NA
bf.norm_approx_logis = NA
bf.norm_approx_beta = NA


null_prob = 0.5
y_prior = 828 # information from a previous study. To get uninformed prior: y_prior = 6.5 and N_prior = 13 will result in rscale = 0.5, y_prior = 0.5 and N_prior = 1 will result in rscale = 1 in the proportionBF function of the BayesFactor package
N_prior = 1560 # information from a previous study. To get uninformed prior: y_prior = 6.5 and N_prior = 13 will result in rscale = 0.5, y_prior = 0.5 and N_prior = 1 will result in rscale = 1 in the proportionBF function of the BayesFactor package
interval = c(0.5,1)

# the below is needed to calculate the rscale parameter for the proportionBF within BayesFactor package
alpha = y_prior + 1
beta = N_prior - y_prior + 1
mean = alpha/(alpha + beta)
var = (alpha*beta)/(((alpha + beta)^2)*(1 + alpha + beta))
prior_mean_lo = qlogis(mean)
prior_sd_lo = sqrt( var / dlogis( prior_mean_lo )^2 )

for(i in 1:length(multip)){
  y = ys[i]
  N = Ns[i]
  
  bf.beta[i] = 1/BF01_beta(y = y, N = N, y_prior = y_prior, N_prior = N_prior, interval = interval, null_prob = null_prob)
  bf.logis[i] = 1/BF01_logis(y = y, N = N, y_prior = y_prior, N_prior = N_prior, interval = interval, null_prob = null_prob)
  bf.norm_approx_logis[i] = 1/BF01_norm_approx_logis(y = y, N = N, y_prior = y_prior, N_prior = N_prior, interval = interval, null_prob = null_prob)
  bf.norm_approx_beta[i] = 1/BF01_norm_approx_beta(y = y, N = N, y_prior = y_prior, N_prior = N_prior, interval = interval, null_prob = null_prob)

}  

par( mfrow = c( 1, 2 ) )

# plot for small sample sizes
plot(ys[1:20],bf.beta[1:20],ty='n',log="y",ylab="Bayes factor",xlab="y",las=1, ylim = c(min(c(bf.beta[1:20], bf.logis[1:20], bf.norm_approx_logis[1:20], bf.norm_approx_beta[1:20])),max(c(  bf.beta[1:20], bf.logis[1:20], bf.norm_approx_logis[1:20], bf.norm_approx_beta[1:20]))))
lines(ys[1:20],bf.beta[1:20],col="black") # BF from the beta distribution
lines(ys[1:20],bf.logis[1:20],col="blue") # normal approcimation of the beta distribution, should match BF from beta distribution if y and N and prior y and N are high enough
lines(ys[1:20],bf.norm_approx_beta[1:20],col="green") # BF from the logistic distribution, same as used in the BayesFactor package in the proportionBF function
lines(ys[1:20],bf.norm_approx_logis[1:20],col="red") # normal approximation of the logistic prior used in the BayesFactor package in the proportionBF function, should match the BFs from logistic distributions closely, but it is not a perfect match

# plot for large sample sizes
plot(ys[290:300],bf.beta[290:300],ty='n',log="y",ylab="Bayes factor",xlab="y",las=1, ylim = c(min(c(bf.beta[290:300], bf.logis[290:300], bf.norm_approx_logis[290:300], bf.norm_approx_beta[290:300])),max(c(  bf.beta[290:300], bf.logis[290:300], bf.norm_approx_logis[290:300], bf.norm_approx_beta[290:300]))))
lines(ys[290:300],bf.beta[290:300],col="black") # BF from the beta distribution
lines(ys[290:300],bf.logis[290:300],col="blue") # normal approcimation of the beta distribution, should match BF from beta distribution if y and N and prior y and N are high enough
lines(ys[290:300],bf.norm_approx_beta[290:300],col="green") # BF from the logistic distribution, same as used in the BayesFactor package in the proportionBF function
lines(ys[290:300],bf.norm_approx_logis[290:300],col="red") # normal approximation of the logistic prior used in the BayesFactor package in the proportionBF function, should match the BFs from logistic distributions closely, but it is not a perfect match







##### demonstration with null effect in the new study and prior distribution based on results of a previous study finding a small positive effect

observed_p = 0.50

multip = 1:300 # to simulate different sample sizes with the same observed p
ys <- 100*observed_p*multip
Ns <- 100*multip

bf.beta = NA
bf.logis = NA
bf.norm_approx_logis = NA
bf.norm_approx_beta = NA

null_prob = 0.5
y_prior = 828 # information from a previous study. To get uninformed prior: y_prior = 6.5 and N_prior = 13 will result in rscale = 0.5, y_prior = 0.5 and N_prior = 1 will result in rscale = 1 in the proportionBF function of the BayesFactor package
N_prior = 1560 # information from a previous study. To get uninformed prior: y_prior = 6.5 and N_prior = 13 will result in rscale = 0.5, y_prior = 0.5 and N_prior = 1 will result in rscale = 1 in the proportionBF function of the BayesFactor package
interval = c(0.5,1)

# the below is needed to calculate the rscale parameter for the proportionBF within BayesFactor package
alpha = y_prior + 1
beta = N_prior - y_prior + 1
mean = alpha/(alpha + beta)
var = (alpha*beta)/(((alpha + beta)^2)*(1 + alpha + beta))
prior_mean_lo = qlogis(mean)
prior_sd_lo = sqrt( var / dlogis( prior_mean_lo )^2 )

for(i in 1:length(multip)){
  y = ys[i]
  N = Ns[i]
  
  bf.beta[i] = 1/BF01_beta(y = y, N = N, y_prior = y_prior, N_prior = N_prior, interval = interval, null_prob = null_prob)
  bf.logis[i] = 1/BF01_logis(y = y, N = N, y_prior = y_prior, N_prior = N_prior, interval = interval, null_prob = null_prob)
  bf.norm_approx_logis[i] = 1/BF01_norm_approx_logis(y = y, N = N, y_prior = y_prior, N_prior = N_prior, interval = interval, null_prob = null_prob)
  bf.norm_approx_beta[i] = 1/BF01_norm_approx_beta(y = y, N = N, y_prior = y_prior, N_prior = N_prior, interval = interval, null_prob = null_prob)

}  

par( mfrow = c( 1, 2 ) )

# plot for small sample sizes
plot(ys[1:20],bf.beta[1:20],ty='n',log="y",ylab="Bayes factor",xlab="y",las=1, ylim = c(min(c(bf.beta[1:20], bf.logis[1:20], bf.norm_approx_logis[1:20], bf.norm_approx_beta[1:20])),max(c(  bf.beta[1:20], bf.logis[1:20], bf.norm_approx_logis[1:20], bf.norm_approx_beta[1:20]))))
lines(ys[1:20],bf.beta[1:20],col="black") # BF from the beta distribution
lines(ys[1:20],bf.logis[1:20],col="blue") # normal approcimation of the beta distribution, should match BF from beta distribution if y and N and prior y and N are high enough
lines(ys[1:20],bf.norm_approx_beta[1:20],col="green") # BF from the logistic distribution, same as used in the BayesFactor package in the proportionBF function
lines(ys[1:20],bf.norm_approx_logis[1:20],col="red") # normal approximation of the logistic prior used in the BayesFactor package in the proportionBF function, should match the BFs from logistic distributions closely, but it is not a perfect match

# plot for large sample sizes
plot(ys[290:300],bf.beta[290:300],ty='n',log="y",ylab="Bayes factor",xlab="y",las=1, ylim = c(min(c(bf.beta[290:300], bf.logis[290:300], bf.norm_approx_logis[290:300], bf.norm_approx_beta[290:300])),max(c(  bf.beta[290:300], bf.logis[290:300], bf.norm_approx_logis[290:300], bf.norm_approx_beta[290:300]))))
lines(ys[290:300],bf.beta[290:300],col="black") # BF from the beta distribution
lines(ys[290:300],bf.logis[290:300],col="blue") # normal approcimation of the beta distribution, should match BF from beta distribution if y and N and prior y and N are high enough
lines(ys[290:300],bf.norm_approx_beta[290:300],col="green") # BF from the logistic distribution, same as used in the BayesFactor package in the proportionBF function
lines(ys[290:300],bf.norm_approx_logis[290:300],col="red") # normal approximation of the logistic prior used in the BayesFactor package in the proportionBF function, should match the BFs from logistic distributions closely, but it is not a perfect match




