# Question 1c)

N = 1000
p = runif(N)
theta = p^(1/3)
mean(theta)
var(theta)

### Question 2

# Part B, C and D

# Pages for 10 weeks printed by 3 persons
y1 = c(86, 74, 79, 70, 84, 65, 57, 49, 70, 79) 
y2 = c(83, 68, 59, 70, 58, 66, 64, 70, 69, 75)
y3 = c(72, 84, 72, 89, 88, 91, 74, 94, 77, 76)

# Used for the mean of y
y.bar = (mean(y1)+mean(y2)+mean(y3))/3 

# Sufficient statistics for likelihoods
y.sum = numeric(3)
y.sum[1] = sum(y1);	
y.sum[2] = sum(y2); 
y.sum[3] = sum(y3)
n = length(y1)

N = 5000
N1 = N+1
B = 100
lambda = matrix(0,ncol=3,nrow=N1)

# Initialize lambda mean based data
lambda[1,] = c(y.bar,y.bar,y.bar) 

# Sample lambda ~ Gamma(y.sum[j]+1, n + 1/y.bar) from  full conditional for  lambda
for(i in 1:N) {
  for(j in 1:3)   lambda[i+1,j] = rgamma(1,y.sum[j]+1, n + 1/y.bar)
}
lambda.b = lambda[B:N,]
mean(lambda.b)
var(lambda.b)

# Separate means by column (for "by row", 2nd arg is 1)
apply(lambda.b,2,mean) 
apply(lambda.b,2,var)

# Compare to MLE's of frequentist. Here are the means from Frequentist statistics for each Yi.
y.sum/n               

# No alpha to plot
# The trace data for each plot shows a good mixture of points and looks very random.
# The histograms show a Gamma distribution.
# The autocorrelation plots go down to 0 very quickly.
# All plots are a good sign of convergence.
par(mfrow=c(3,3))     
for(j in 1:3) plot(lambda[,j],pch=".",ylab=paste("lambda[",j,"]"))
for(j in 1:3) hist(lambda.b[,j],pch=".",main=paste("Histogram of theta[",j,"]"))
for(j in 1:3) acf(lambda[,j],main=paste("Series theta[",j,"]"))

# Credible intervals for each lamba j
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

### Question 3

# Log of full conditional for alpha
flog.alpha = function(alpha,lambda) {     
  f = -3*lgamma(alpha)+alpha*log(lambda[1]*lambda[2]*lambda[3])-alpha/100
  return(f)
}

# Pages for 10 weeks
y1 = c(86, 74, 79, 70, 84, 65, 57, 49, 70, 79) 
y2 = c(83, 68, 59, 70, 58, 66, 64, 70, 69, 75)
y3 = c(72, 84, 72, 89, 88, 91, 74, 94, 77, 76)

# Used for mean of y
y.bar = (mean(y1)+mean(y2)+mean(y3))/3

# Sufficient statistics for likelihoods
y.sum = numeric(3)
y.sum[1] = sum(y1);	y.sum[2] = sum(y2); y.sum[3] = sum(y3)
n = length(y1)

N = 5000
N1 = N+1
B = 1000	
alpha = c()
lambda = matrix(0,ncol=3,nrow=N1)
alpha[1] = y.bar          

# Initialize lambda mean based on alpha prior
lambda[1,] = c(y.bar,y.bar,y.bar)               
accept.alpha = 0
d = 10

# Sample lambda ~ Gamma(alpha[i]+y.sum.j, n + 1), may need to tune sd 
for(i in 1:N) {
  for(j in 1:3)   lambda[i+1,j] = rgamma(1,alpha[i]+y.sum[j], n + 1)
  alpha.star = -1
  while(alpha.star < 0) { 
    alpha.star = rnorm(1,alpha[i],d)
  }
  f.star = flog.alpha(alpha.star,lambda[i+1,])
  fi = flog.alpha(alpha[i],lambda[i+1,])
  alpha.prob = min(1,exp(f.star - fi))
  if(runif(1) < alpha.prob) {
    alpha[i+1] = alpha.star
    accept.alpha = accept.alpha + 1
  }else alpha[i+1] = alpha[i]
}
accept.alpha/N
alpha.b = alpha[B:N]
lambda.b = lambda[B:N,]
mean(alpha.b)
var(alpha.b)
mean(lambda.b)
var(lambda.b)

# Separate means by column (for "by row", 2nd arg is 1)
apply(lambda.b,2,mean) 
apply(lambda.b,2,var) 

# Compare to MLE's
y.sum/n               
par(mfrow=c(3,4))
plot(alpha,pch=".")
for(j in 1:3) plot(lambda[,j],pch=".",ylab=paste("lambda[",j,"]"))
hist(alpha.b)
for(j in 1:3) hist(lambda.b[,j],pch=".",main=paste("Histogram of theta[",j,"]"))
acf(alpha)
for(j in 1:3) acf(lambda[,j],main=paste("Series theta[",j,"]"))

# Credible intervals
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

