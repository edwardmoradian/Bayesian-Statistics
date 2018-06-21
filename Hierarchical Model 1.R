### Hierarchical model
### The data is derived from an exponential distribution of unknowwn lambda
### data: 5.41  3.95  4.00  2.46  3.07 21.24  3.49  2.67  2.29 16.98
### 
### It is thought that this data comes from a process (1/5)exp(-lambda/5)
### The following proper priors reflect that intuition:
### A prior is put on lambda: 
###    lambda ~ lambda^(alpha.0-1)*exp(-beta*lambda)
###    lambda ~ gamma(alpha.0,beta) 
###    where alpha.0 = .2
### A prior is put on beta: beta ~ alpha*exp(-alpha*beta) ~ Exp(alpha)
### A prior is put on alpha: alpha ~ exp(-alpha) ~ Exp(1)
###
### The hierarchy is
###        X[i] ~ Exp(lambda)
###        lambda ~ Gamma(alpha.0,beta)
###        beta ~ Exp(alpha)
###        alpha ~ Exp(1)
###
### Work out the full conditionals for lambda, theta, alpha
### to from 
###  
### lamda ~ Gamma(10+alpha.0,sum(x)+beta) 
### beta ~ Gamma(alpha.0+1,lambda + alpha)
### alpha ~ Gamma(2,beta+1)
###
### They should be recognizable distributions
###
### Do a Gibbs sampling to find the means and standard deviations of lambda, beta, alpha
###
# Initialize

x = c(5.41,  3.95,  4.00,  2.46,  3.07,  21.24,  3.49,  2.67,  2.29, 16.98)
alpha.0 = .2
x.bar = mean(x)
sumx = sum(x)
n = length(x)
lambda = c()
lambda[1] = 1/x.bar
lambda[1]
beta = c()
beta[1] = 1
alpha = c()
alpha[1] = 1
N = 10000
Nb = 2000
N1 = N+1

for(i in 1:N) 
{
	# fill in the Gibbs sampling for each parameter
	lambda[i+1] = rgamma(1,10 + alpha.0,sum(x)+beta[i])
	beta[i+1] = rgamma(1,alpha.0+1,lambda[i+1]+alpha[i])
  alpha[i+1] = rgamma(1,2,beta[i+1]+1)
}

mean(lambda[Nb:N])
sd(lambda[Nb:N])
mean(beta[Nb:N])
sd(beta[Nb:N])
mean(alpha[Nb:N])
sd(alpha[Nb:N])
par(mfrow=c(3,3))
hist(lambda[Nb:N])
hist(beta[Nb:N])
hist(alpha[Nb:N])
plot(lambda,type='p',pch='.')
plot(beta,type='p',pch='.')
plot(alpha,type='p',pch='.')
acf(lambda)
acf(beta)
acf(alpha)

### Now suppose we put a flat prior on lambda: f(lambda) = 1 
###  Do the Gibbs sampling for lambda

lambda = c()
lambda[1] = 1/x.bar

for(i in 1:N) 
{
	lambda[i+1] = rgamma(1,11,sum(x))
}

mean(lambda[Nb:N])
sd(lambda[Nb:N])
par(mfrow=c(3,1))
hist(lambda)
plot(lambda,type='p',pch='.')
acf(lambda)
	
### What is the difference between the two means for lambda?
### Can the priors be justified? Are they informative?

# The priors are uninformative.  They seem to reduce the mean of lambda
# more than expected. From 0.168 to 0.152.
# The intention was to give students work on a hierarchical model.
# 



