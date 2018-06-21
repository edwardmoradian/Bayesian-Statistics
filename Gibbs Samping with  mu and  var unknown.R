# data for example
# from y = rnorm(20,15,2)
y = c(11.54053,16.60317,15.99318,15.38895,13.18257,18.03025,15.61183,14.42232,
      16.32630,14.43119,11.34286,16.12436,13.32910,14.78058,14.88083,15.48779,
      17.00922,14.08266,18.12429,13.31570)
n = length(y)
# intialize paramters 
mu = c()            # make mu a vector
y.bar = mean(y)
mu[1] = mean(y)     # initialize mu
sigma.2 = c()       # make sigma.2 a vector
sigma.2[1] = var(y) # initialize sigma.2
# set run parameters
N = 5000
N1 = N + 1
B = 1000
for(i in 1:N1) {
	mu[i+1] = rnorm(1,y.bar,sqrt(sigma.2[i]/n))
	sigma.2[i+1] = 1/rgamma(1,n/2,1/2*sum((y-mu[i+1])^2)) 
}
#    compute means and variances excluding burn in iterations
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
# check with frequentist values
mean(y)
var(y)/n
t.alpha.2 =qt(.975,19)
U = mean(y) + t.alpha.2*sqrt(var(y)/n)
L = mean(y) - t.alpha.2*sqrt(var(y)/n)
msg = paste("(",L,",",U,")")
print(msg,quote=FALSE)
print(c("sigma.2.freq = ",var(y)),quote=FALSE)



