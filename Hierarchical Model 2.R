flog.alpha = function(alpha,beta,lambda) {
	f = 3*alpha*log(beta)+624*log(alpha)-3*lgamma(alpha)-2.5*alpha + alpha*sum(log(lambda))
return(f)
}
# Note: For our homework 3, the above is changed. 

#cookie data for brands chips, niteowl, safeway
chips = c(28, 25, 22, 28, 25, 22, 21, 27, 23, 22, 24, 22, 23, 25, 22, 27, 
          24, 22, 26, 22, 22, 18, 23, 24, 21, 19, 22, 25, 21, 21, 26)
niteowl = c(28, 25, 22, 28, 25, 30, 34, 27, 23, 31, 24, 22, 23, 25, 22, 27,
            24, 30, 26, 22, 22, 18, 23, 24, 30, 19, 22, 25, 21, 30, 26)
safeway = c(28, 25, 22, 28, 25, 22, 21, 27, 23, 31, 24, 22, 23, 25, 22, 27, 
            24, 22, 26, 22, 22, 18, 23, 24, 30, 19, 22, 25, 21, 21, 26)

y.sum = numeric(3)
y.sum[1] = sum(niteowl) 
y.sum[2] = sum(chips) 
y.sum[3] = sum(safeway) 
n = numeric(3)
n[1] = length(niteowl)
n[2] = length(chips)
n[3] = length(safeway)

# important statistics
y.sum
n

N = 5000
N1 = N+1
B = 100
alpha = c()
lambda = matrix(0,ncol=3,nrow=N1)
alpha[1] = 250                           # prior mean
beta = 10                                # fixed beta for prior 
lambda[1,] = c(25,25,25)                 # initialize lambda mean based og alpha prior, first row of alpha
accept.alpha = 0
flg.int = 1
for(i in 1:N) {
	#sample lambda ~ Gamma(alpha[i]+y.sum.j, 10 + n.j)
	for(j in 1:3)   lambda[i+1,j] = rgamma(1,alpha[i]+y.sum[j], 10+n[j])
	alpha.star = -1
	while(alpha.star < 0) {
  		alpha.star = rnorm(1,alpha[i],sd=15)   # may need to tune sd 
  	}
	f.star = flog.alpha(alpha.star,beta,lambda[i,])
	fi = flog.alpha(alpha[i],beta,lambda[i,])
	alpha.prob = min(1,exp(f.star - fi))
	if(runif(1) < alpha.prob) {
		alpha[i+1] = alpha.star
		accept.alpha = accept.alpha + 1
	}else alpha[i+1] = alpha[i]    #use current iteration as the new iteration
}
accept.alpha/N
alpha.b = alpha[B:N]
lambda.b = lambda[B:N,]
mean(alpha.b)
var(alpha.b)
mean(lambda.b)
var(lambda.b)
apply(lambda.b,2,mean)   #separate means by column (for "by row", 2nd arg is 1)
y.sum/n                #compare to MLE's
par(mfrow=c(3,4))
plot(alpha,pch=".")
for(j in 1:3) plot(lambda[,j],pch=".",ylab=paste("lambda[",j,"]"))
hist(alpha.b)
for(j in 1:3) hist(lambda.b[,j],pch=".",main=paste("Histogram of theta[",j,"]"))
acf(alpha)
for(j in 1:3) acf(lambda[,j],main=paste("Series theta[",j,"]"))
#credible intervals
lam1 = lambda.b[,1]
lam1 = sort(lam1)
nn = length(lam1)
L = as.integer(.025*nn)
U = as.integer(.975*nn)
print(c("CI for lambda 1:",lam1[L],lam1[U]),quote=FALSE)

lam2 = lambda.b[,2]
lam2 = sort(lam2)
print(c("CI for lambda 2:",lam2[L],lam2[U]),quote=FALSE)

lam3 = lambda.b[,3]
lam3 = sort(lam3)
print(c("CI for lambda 3:",lam3[L],lam3[U]),quote=FALSE)


