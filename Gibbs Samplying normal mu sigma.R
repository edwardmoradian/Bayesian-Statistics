# code for gibbs sampling example for normal likelihood with mu and tau unknown.
# initialize data
# 10 x values drawn from rnorm(10,10,1), using sigma = 1
x= c(11.001447,9.249269,9.437698,8.922587,11.971188,10.519513,8.699515,10.878161,11.350819,9.834667)
n = length(x)         
s.0 = 0.5      # set s.0 for prior on mu
m = 9.5        # set m for prior on mu
alpha = 3      # set alpha for prior on sigma.2
beta = 2       # set beta for prior on sigma.2
mu = c()       # make mu a vector
mu[1] = m      # intialize mu  
sigma.2 = c()  # make sigm.2 a vector
sigma.2[1] = beta/(alpha-1) # intialize sigma.2
# compute constants needed for iterations loop
x.bar = mean(x) # mean of x
mu.cons = (n*x.bar+s.0*m)/(n+s.0) # mu for full conditional 
alpha.sigma.2 = n/2+ alpha        # alpha for inverse gamma for sigma.2
beta.sigma.2 = (1/2)*sum((x-x.bar)^2)+(1/2)*(n*s.0*(x.bar-m)^2)/(n+s.0)+beta
#               beta for inverse gamma for sigma.2
#set loop counts
N = 10000
B = 1000
N1 = N+1       # iteration count - one above N
#begin iteration loop
for(i in 1:N1) {
	# get next iteration of mu
	mu[i+1] = rnorm(1,mu.cons,sqrt(sigma.2[i]/(n+s.0))) # Gibbs sample mu
      # get next iteration of tau
	sigma.2[i+1] = 1/rgamma(1,alpha.sigma.2,beta.sigma.2) # Gibbs sample sigma.2
}
#     compute mean and variances of mu and sigma.2 exculding burn in iterations
mean(mu[B:N1])
var(mu[B:N1])
mean(sigma.2[B:N1])
var(sigma.2[B:N1])
par(mfrow=c(2,3))    # graph histograms, traces and autocorrelations
hist(mu[B:N1])       # should be t distribution
plot(mu,type='p',pch='.') # shows good mixing
acf(mu)              # shows little correlation of iterates
#
hist(sigma.2[B:N1])      # should be a gamma distribution
plot(sigma.2,type='p',pch='.')   # shows good mixing
acf(sigma.2)             # shows little correlation of iterates
# credible intervals
# credible interval for mu
mu.cr = mu[B:N1]
mu.cr = sort(mu)
sigma.2.cr = sigma.2[B:N1]
sigma.2.cr = sort(sigma.2)
L = as.integer((N1-B)*.025)
U = as.integer((N1-B)*.975)
msg.cred.mu = paste(c("(",mu.cr[L],",",mu.cr[U],")"))
print(msg.cred.mu,quote=FALSE)
#credible interval for tau
msg.cred.sigma.2 = paste(c("(",sigma.2.cr[L],",",sigma.2.cr[U],")"))
print(msg.cred.sigma.2,quote=FALSE)
#check with data
#frequentist values
x.bar         # check with mean of mu, expect to be a little greater
var(x)/n      # check with the variance of mu, expect to be similar
t.alpha.2 = qt(.975,9)
L.frq = x.bar-t.alpha.2*sd(x)/sqrt(n)
U.frq = x.bar+t.alpha.2*sd(x)/sqrt(n)
# confidence interval
msg.freq.mu = paste(c("(",L.frq,",",U.frq,")"))
print(msg.freq.mu,quote=FALSE)
# further checks for sigma.2 mean, sigma.2 variance
sigma.2.data = var(x)
sigma.2.data      # check with mean of sigma.2, expect to be similar
beta.sigma.2/(alpha.sigma.2-1)   # check with mean of marginal distribution of sigma.2,
#                                 expect to be close
beta.sigma.2^2/((alpha.sigma.2-1)^2*(alpha.sigma.2-2)) # check with variance of the marginal distribution of
#                                                        sigma.2, expect to be close



