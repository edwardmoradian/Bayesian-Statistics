### Comparing proportion tests for a coin flip data.
### A) First use draws from the distribution simulating
### iterations from a Gibbs sampling or M-H sampling
###
### B) Next use the distribution itself.
### 
### Consider a coin flip experiment.  
### The coin is flipped 100 times to get 40 heads.  
### Assume a flat prior of Beta(1,1).
###
### Determine the Beta distribution.
###
### Beta(?,?)  <-fill alpha and beta
###
### Make N = 10000 draws of theta from this distribution to 
### simulate Gibbs sampling. Call the variable theta.
###
N = 10000
alpha = 41
beta = 61
theta = rbeta(N,alpha,beta) 

## Use the draws to find the mean of the distribution.
## Let the mean be theta.bar 
###
theta.bar = mean(theta)
theta.bar
###
### Make a histogram hist(theta,xlim=c(0,1) and
### draw a vertical line for the mean on the histogram 
### using abline(v=theta.bar,col='red)
###
par(mfrow=c(1,1))
hist(theta,xlim=c(0,1))
abline(v=theta.bar,col='red')
###
### Do the hypothesis test: Use the iterations
###     H0: theta => .5  vs H1: theta < .5
###     at significance level alpha = .05
### That is, determine the probability the distribution is greater than
### .5 by finding the probability P(theta < .5) using the draws
### from the Beta distribution.  [Hint: ct = sum(theta >= .5)]
###  Name this probability P.gt.half.
ct = sum(theta >= .5)                # ct is count of theta > .5
P.gt.half = ct/N      # probability greater than .5
P.gt.half
if( P.gt.half <= .05) {
	print("reject H0",quote=FALSE)
}else print("do not reject H0",quote=FALSE)
###
### Do the hypothesis test: Use the Distribution
###     H0: theta = .5  vs H1: theta != .5
###     at significance level alpha = .05
### That is, determine the 95% CI for draws for theta and see if it 
### includes .5. Name upper and lower values theta.L and theta,U
###
theta = sort(theta)
L = as.integer(.025*N)
U = as.integer(.975*N)
theta.L = theta[L]
theta.U = theta[U]
theta.L
theta.U
if( theta.L <= .5 && theta.U >= .5 ) {
	print("do not reject H0",quote=FALSE)
}else print("reject H0",quote=FALSE)
###
### Now use the probaility distribution Beta(alpha,beta) to do these test.
###
### Do the hypothesis test:
###     H0: theta <= .5  vs H1: theta > .5
###     at significance level alpha = .05
### That is, determine the probability the distribution is greater than
### .5 by using pbeta to get the probabiity greater than .5.
###
P.gt.half = 1 - pbeta(.5,alpha,beta)  # accumulates up .5
P.gt.half
if( P.gt.half <= .05) {
	print("reject H0",quote=FALSE)
}else print("do not reject H0",quote=FALSE)
###
### Do the hypothesis test:
###     H0: theta = .5  vs H1: theta != .5
###     at significance level alpha = .05
### That is, determine the 95% CI for distribution theta using qbeta
### and see if it includes .5. 
###
theta.L = qbeta(.025,alpha,beta)
theta.U = qbeta(.975,alpha,beta)
theta.L
theta.U
if( theta.L <= .5 && theta.U >= .5 ) {
	print("do not reject H0",quote=FALSE)
}else print("reject H0",quote=FALSE)
###
### How do the probabilities compare for the first test?
### How do the intervals compare for the second test?
