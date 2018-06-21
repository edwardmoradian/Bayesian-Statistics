### 1 parameter problem

f=function(theta) 
{ 
  dbeta(theta,7,5)
}

# loop values
N = 5000
N1 = N + 1
B = 1000
accept.ct = 0   # acceptance count
theta = c()     # make theta a vector
theta[1] = runif(1)  # intialize theta - we don't know the function

for(i in 1:N1) 
{
	theta.star = 2
while (theta.star > 1 || theta.star < 0 ) 
{   
	theta.star = runif(1,theta[i]-.2,theta[i]+.2) # symmetric proposal uniform, draw from uniform
}
	f.star = f(theta.star)   # compute beta function for theta.star
	fi = f(theta[i])         # compute beta function for theta[i]
	alpha = min(1,f.star/fi) # alpha is value compared to runif(1), calculate the acceptance probability
if(runif(1) <= alpha ) 
{
	theta[i+1] = theta.star   # Accept the theta value
	accept.ct = accept.ct + 1
}
else 
{
  theta[i+1] = theta[i] # When not accepting, set next value to the current value
}
}

# acceptance count - helps to set proposal function
accept.ct/N 

# compute mean and variance ecluding burn in iterations

mean(theta[B:N])
var(theta[B:N])

# credible interval

theta.cr = theta[B:N]
theta.cr = sort(theta.cr)
U = as.integer(.975*(N-B))
L = as.integer(.025*(N-B))
print(c("(",round(theta.cr[L],4),",",round(theta.cr[U],4),")"),quote=FALSE)
par(mfrow=c(3,1))
hist(theta[B:N])
plot(theta,type='p',pch='.')
acf(theta)

# check actual values

theta.mean.actual = round(7/(7+5),4)
var.actual = round(5*7/(12^2*13),4)
print(c("Actual theta mean and variance =",
	theta.mean.actual,var.actual),quote=FALSE)


