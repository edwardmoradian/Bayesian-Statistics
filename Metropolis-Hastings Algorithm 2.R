### 2 parameter problem

logf = function(mu,sigma.2,y,n)  #log of the posterior distribution - values may be small or large
{
  logf = -(n/2+1)*log(sigma.2)-1/(2*sigma.2)*sum((y-mu)^2)
	return(logf)
}

# data for example from y = rnorm(20,15,2)
y = c(11.54053,16.60317,15.99318,15.38895,13.18257,18.03025,15.61183,14.42232,
      16.32630,14.43119,11.34286,16.12436,13.32910,14.78058,14.88083,15.48779,
      17.00922,14.08266,18.12429,13.31570)
n = length(y)

# intialize parameters 
mu = c()            # make mu a vector
mu[1] = mean(y)     # initialize mu
sigma.2 = c()       # make sigma.2 a vector
sigma.2[1] = var(y) # initialize sigma.2

# set run parameters
N = 5000
N1 = N + 1
B = 1000
accept.mu.ct = 0
accept.sigma.ct = 0
d = sd(y)/sqrt(n)  # select a standard deviation for drawing mu.star from proposal distribution

for(i in 1:N1) 
{
	# M-H on mu
	mu.star = rnorm(1,mu[i],d) # draw from normal distribution (symmetric proposal distribution)
	logf.star = logf(mu.star,sigma.2[i],y,n)
	logfi = logf(mu[i],sigma.2[i],y,n)
	alpha = min(1,exp(logf.star-logfi))   # for log of posterior, compare using exp()

if(runif(1) <= alpha) 
{
		mu[i+1] = mu.star
		accept.mu.ct = accept.mu.ct + 1
}
else mu[i+1] = mu[i]
{
	# M-H on sigma.2
	sigma.2.star = runif(1,sigma.2[i]/2,2*sigma.2[i])    # draw from asymmetric proposal distribution
	logf.star = logf(mu[i+1],sigma.2.star,y,n)
	logfi = logf(mu[i+1],sigma.2[i],y,n)
	logf.star.q = logf.star-(log(2)-log(3*sigma.2[i])) # divide by proposal probability
	logfi.q = logfi-(log(2)-log(3*sigma.2.star))       # divide by proposal probability
	alpha = min(1,exp(logf.star.q-logfi.q))          
}
if(runif(1) <= alpha) 
{
		sigma.2[i+1] = sigma.2.star
		accept.sigma.ct = accept.sigma.ct + 1
}
else sigma.2[i+1] = sigma.2[i]
}

# acceptance probabilities
accept.mu.ct/N
accept.sigma.ct/N

# compute means and variances excluding burn in iterations
mean(mu[B:N])
var(mu[B:N])
mean(sigma.2[B:N])
var(sigma.2[B:N])

# credible interval for mu
mu.cr = mu[B:N]
mu.cr = sort(mu.cr)
U = as.integer(.975*(N-B))
L = as.integer(.025*(N-B))
print(c("(",round(mu.cr[L],4),",",round(mu.cr[U],4),")"),quote=FALSE)

# credible interval for sigma.2
sigma.2.cr = sigma.2[B:N]
sigma.2.cr = sort(sigma.2.cr)
U = as.integer(.975*(N-B))
L = as.integer(.025*(N-B))
msg = paste("(",round(sigma.2.cr[L],4),",",round(sigma.2.cr[U],4),")")
print(msg,quote=FALSE)

#  histograms, traces and autocorrelations
par(mfrow=c(2,3))
hist(mu[B:N],freq=FALSE)
plot(mu,type='p',pch='.')
acf(mu)
hist(sigma.2[B:N],freq=FALSE)
plot(sigma.2,type='p',pch='.')
acf(sigma.2)

#dev.off()
# check with frequentist values
mean(y)
var(y)/n
t.alpha.2 =qt(.975,19)
U = mean(y) + t.alpha.2*sqrt(var(y)/n)
L = mean(y) - t.alpha.2*sqrt(var(y)/n)
msg = paste("(",L,",",U,")")
print(msg,quote=FALSE)




