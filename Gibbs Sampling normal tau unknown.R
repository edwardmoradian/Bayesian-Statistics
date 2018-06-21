# code for gibbs sampling example for normal likelihood 
# with mu known and tau unknown.
# initialize data
# 10 x values drawn from rnorm(10,10,1), using sigma = sqrt(1/tau) for normal functions ==> tau = 1
x= c(11.001447,9.249269,9.437698,8.922587,11.971188,10.519513,8.699515,10.878161,11.350819,9.834667)
n = length(x)  
x.bar = mean(x)
alpha = 1
b = .5  # constant for prior of beta --> b*exp(-b*beta)
beta = c()  # make beta a vector
beta[1] = 1 # initialize hyperparameter beta
tau = c()   # make tau a vector
mu = x.bar  # set mu to x.bar, it is constant for this application
tau[1] = alpha/beta[1] # intialize tau
# compute constants needed for iteration loop
alpha.tau = n/2+ alpha    # alpha for gamma for tau
beta.tau.0 = (1/2)*sum((x-mu)^2)   # beta for gamma for tau
alpha.beta = alpha + 1    # alpha for gamma for beta 
#set loop counts
N = 10000
B = 1000
N1 = N+1
#begin iteration loop
for(i in 1:N1) {
     	tau[i+1] = rgamma(1,alpha.tau,beta.tau.0 + beta[i])  # get next iteration of tau
      beta[i+1] = rgamma(1,alpha.beta,b + tau[i+1])  # get next iteration of beta
}
#   compute mean and variances excluding burn in iterations 
mean(tau[B:N1])
var(tau[B:N1])
mean(beta[B:N1])
var(beta[B:N1])
par(mfrow=c(2,3))
hist(tau[B:N1])       # gamma distribution
plot(tau,type='p',pch='.') # shows good mixing
acf(tau)              # shows little correlation of iterates
#
hist(beta[B:N1],breaks=20)      # should be a gamma distribution
plot(beta,type='p',pch='.')   # shows good mixing
acf(beta)             # shows little correlation of iterates
# credible intervals
# credible interval for tau
tau.cr = tau[B:N1]
tau.cr = sort(tau.cr)
beta.cr = beta[B:N1]
beta.cr = sort(beta.cr)
L = as.integer((N1-B)*.025)
U = as.integer((N1-B)*.975)
msg.cred.tau = paste(c("(",tau.cr[L],",",tau.cr[U],")"))
print(msg.cred.tau,quote=FALSE)
#credible interval for beta
msg.cred.beta = paste(c("(",beta.cr[L],",",beta.cr[U],")"))
print(msg.cred.beta,quote=FALSE)
#check with data
#frequentist point estimate
tau.pt.estimate.freq = (n-1)/sum((x-x.bar)^2)
tau.pt.estimate.freq




