### Gibbs and Metropolis-hastings sampling for a normal distribution.
###
### Assignment - 1) add the code for the Gibbs and M-H sampling.
###              2) run the code to get the means and standard deviations for mu, sigma.2,
###  	            graphs of histograms, traces, autocorrelations, and CI's
###              3) experiment with d, the proposal distribution 
###                 standard deviation to get the best mixing
###              
# Data: 16.22 10.70 16.23 16.63 10.96 15.18 16.67 17.33 18.72 16.12
#
# f(mu|sigma.2,x) ~ N(x.bar,sigma.2/n)
# f(sigma.2|mu,x) proporional to sigma.2^(-n/2-1)*exp[-1/(2*sigma.2)*(sum((mu-x)^2)]
# Where f(mu) = 1 and f(sigma.2) = 1/sigma.2
#
f.sigma.2 = function(sigma.2,mu,n,x) 
{
	f = sigma.2^(-(n/2+1))*exp(-1/(2*sigma.2)*sum((mu-x)^2))
	return(f)
}
x = c(16.22, 10.70, 16.23, 16.63, 10.96, 15.18, 16.67, 17.33, 18.72, 16.12)
n = length(x)
x.bar = mean(x)
mu = c()
mu[1] = x.bar
sigma.2 = c()
sigma.2[1] = var(x)
N = 10000
Nb = 2000
N1 = N + 1
accept = 0     # counter for acceptance ratio
d = 5 # standard deviation - experiment with this value for best mixing
for(i in 1:N) 
{
  mu[i+1] = rnorm(1,x.bar,sqrt(sigma.2[i]/n)) # Gibbs sample mu[i+1] from Normal distribution
	sigma.2.star = -1 # Use normal proposal distribution to draw sigma.2.star
while(sigma.2.star < 0) 
{
  sigma.2.star=rnorm(1,sigma.2[i],d) # Remember make sure it is > 0.  A while loop can be helpful here.
} 
	f.star= f.sigma.2(sigma.2.star,mu[i+1],n,x) # Call function f.sigma.2 to compute full conditional for sigma.2.star
	f.i = f.sigma.2(sigma.2[i],mu[i+1],n,x) # Call function f.sigma.2 to compute full conditional for sigma.2[i]
	min(1,f.star/f.i) #compute the acceptance probability
if( runif(1) <= min(1,f.star/f.i) )  # Use if-else statement for comparison to runif(1) Uniform distribution
{ 
  sigma.2[i+1] = sigma.2.star  # to set sigma.2[i+1] to sigma.2.star else sigma.2[i+1] to sigma.2[i]
  accept = accept+1   # Increment acceptance counter as required 
}
	else sigma.2[i+1]=sigma.2[i]  # When if statement is a not executed, else set sigma.2[i+1=sigma.2[i]
}

accept/N            # acceptance ratio
mean(mu[Nb:N])      # mu and sigma.2 means and standard deviations
sd(mu[Nb:N])
mean(sigma.2[Nb:N])
sd(sigma.2[Nb:N])
par(mfrow=c(3,2))   # graphs of histograms, traces, and autocorrelations
hist(mu)
hist(sigma.2)
plot(mu,type='p',pch='.')
plot(sigma.2,type='p',pch='.')
acf(mu)
acf(sigma.2)
##################### credible intervals
mu.b = mu[Nb:N]
mu.b = sort(mu.b)
sigma.2.b = sigma.2[Nb:N]
sigma.2.b = sort(sigma.2.b)
Ns = length(mu.b)
L = as.integer(.025*Ns)
U = as.integer(.975*Ns)
msg = paste("CI for mu: ",round(mu.b[L],3),",",round(mu.b[U],3))
print(msg,quote=FALSE)
msg = paste("CI for sigma.2: ",round(sigma.2.b[L],3),",",round(sigma.2.b[U],3))
print(msg,quote=FALSE)

### Note: You could have used Gibbs sampling for sigma.2.
###       The inverse gamma distribution associated with sigma.2 is 
###          Inv_gamma(alpha,beta) where alpha = n/2, beta = 1/2*sum((mu[i+1]-x)^2)
###       To sample use sigma.2[i+1] = 1/rgamma(alpha,beta)
###       Try this on your own time.

		
