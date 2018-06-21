# Hierarchical Example 2. 
va1 = c(236,76,55,59,85,12,71,35,122,131,103,22,8,120,104,224,369,69,5,28,306,252,69,95,66,96,100,53,147,99,     
         58,100,145,9,45,67,13,116,125,119,154,117,61,257,137,106,93,273,205,164,115,333,141,26,130,78,227,261,8,174,138,31,48,109,52,100,201,7,49,257,
         21,482,50,77,115,71,245,106,165,44,48,178,171,135,65,320,19,283,210,224,121,86,19,75,108,24,2,65,149,16,212,30,109,127,164,8,135,69,177,108,
        116,97,136,117,113,419,102,187,314,343,18,88,115,37,55,188,91,68,44,138,296,329,207,172,107,159,128,117,96,125,21,80,260,725,193,104,254,150,75,272,
        637,323,48,26,141,163,359,115,113)    
va2 = c(520,186,134,162,175,65,158,63,350,367,201,32,21,366,274,398,812,187,15,124,651,694,116,338,153,328,147,96,451,283,
        195,309,399,25,175,236,26,283,213,336,516,285,165,827,461,211,250,713,497,562,276,506,412,40,361,165,357,465,16,315,247,46,138,193,95,197,445,13,144,603,
         55,773,105,170,246,83,451,367,390,76,295,256,411,258,165,807,41,540,495,535,261,160,42,132,337,33,8,195,417,17,542,40,654,313,472,21,255,125,427,260,
        438,304,304,221,168,1093,313,364,611,822,43,151,370,89,191,385,288,189,71,383,686,706,392,270,252,437,433,336,193,232,38,211,665,1214,
        651,243,489,549,217,659,1013,710,241,75,383,454,1126,294,279)   

va = cbind(va1,va2,deparse.level=0)    # va1 is Y[i], va2 is N[i]
alpha = numeric(N1)      # make a vector for alpha
beta = numeric(N1)       # make a vector for beta
n = length(va[,1])       # number of va hospitals
va
theta = matrix(0,ncol=n,nrow=N1)        # make matrix for theta
alpha[1] = 10                           # initialize prior alpha
beta[1] = 10                            # initialize prior beta
theta[1,] = c(rep(.5,n))                # intialize theta to 0.5
# Set acceptance counts and loop counts
accept.alpha = 0
accept.beta = 0
N = 5000
N1 = N+1
B = 500
for(i in 1:N) {
	# M-H for alpha
	log.theta.sum = sum(log(theta[i,]))    # log of the product of the theta's
	log.theta.minus.sum = sum(log(1-theta[i,]))  # log of the product of the (1-theta)'s
	alpha.star = rnorm(1,alpha[i],sd=1)   
	if(alpha.star > 0) {
		f.star = n*lgamma(alpha.star+beta[i])-n*lgamma(alpha.star)+alpha.star*log.theta.sum-alpha.star/10
		fi = n*lgamma(alpha[i]+beta[i])-n*lgamma(alpha[i])+alpha[i]*log.theta.sum-alpha[i]/10
		alpha.prob = min(1,exp(f.star-fi)) # for division of full conditionals, subtract logs raised to exp
	      if(runif(1) < alpha.prob ) {
		      alpha[i+1] = alpha.star
		      accept.alpha = accept.alpha + 1
		}
		else alpha[i+1] = alpha[i]
	}
	else alpha[i+1] = alpha[i]
      # M-H for beta
  	beta.star = rnorm(1,beta[i],sd=1)   
  	if(beta.star > 0) {
    		f.star = n*lgamma(alpha[i+1]+beta.star)-n*lgamma(beta.star)+beta.star*log.theta.minus.sum-beta.star/10
        	fi = n*lgamma(alpha[i+1]+beta[i])-n*lgamma(beta[i])+beta[i]*log.theta.minus.sum-beta[i]/10
		beta.prob = min(1,exp(f.star-fi))   # for division of full conditionals, subtract logs raised to exp
	      if(runif(1) < beta.prob) {
      		beta[i+1] = beta.star
      		accept.beta = accept.beta + 1
    		}
    		else beta[i+1] = beta[i]
  	}
  	else beta[i+1] = beta[i]
	# gibbs sample theta
  	#sample theta.j ~ Beta(y.j + alpha[i+1], n.j - y.j + beta[i+1])
  	theta[i+1,] = rbeta(n, va[,1]+alpha[i+1], va[,2]-va[,1]+beta[i+1])
}
accept.alpha/N      # acceptance ratio for alpha
accept.beta/N       # acceptance ratio for beta
alpha.b = alpha[B:N]
beta.b = beta[B:N]
theta.b = theta[B:N,]
par(mfrow=c(3,2))

# histograms for alpha, beta, and first 4 theta
hist(alpha.b)
hist(beta.b)
for(j in 1:4) hist(theta.b[,j],main=paste("Histogram of theta.b[",j,"]"))

#traces for alpha, beta, and first 4 theta
plot(alpha,pch=".")
plot(beta,pch=".")
for(j in 1:4) plot(theta[,j],pch=".",ylab=paste("theta[",j,"]"))

# autocorrelations for alpha, beta, and first 4 theta
acf(alpha)
acf(beta)
for(j in 1:4) acf(theta[,j],main=paste("Series theta[",j,"]"))

mean(alpha.b)
var(alpha.b)
mean(beta.b)
var(beta.b)
mean(theta.b)     # overall mean for all hospitals
prior.mean = mean(alpha.b)/(mean(alpha.b)+mean(beta.b))
prior.mean

round(apply(theta.b,2,mean),3)     #separate theta means for each hospital (2 = column, mean = function)
round(va[,1]/va[,2],3)             #compare to MLE's given as Y[i]/N[i]
# the theta means are shrunk towards the a combination 
# of the data mean (approx =.434) and the prior mean (approx = .434)
par(mfrow=c(1,1))
plot(va[,1]/va[,2],apply(theta.b,2,mean),xlab="MLE",ylab="posterior mean") # (x,y)=(MLE[i],theta[i]) ith hopital
abline(0,1)   # line y=x

plot(va[,1]/va[,2],apply(theta.b,2,mean),col=1+(va[,2]<100))  # same graph, hospital with N[i]<100 in red
abline(0,1)               #color 1 is black, color 2 is red
abline(h=mean(theta.b),col="green")    # horizontal line at y = overall theta mean

#probability hospital 1 had a higher dissatisfaction rate than hospital 2
mean(theta.b[,1]>theta.b[,2])

#probability hospital 1 rates in the top half of dissatisfaction rates
theta.rank = apply(theta.b,1,rank)   #compute rank within rows (iterations)
dim(theta.rank)    # rows and columns
mean(theta.rank[1,])  # mean rank of hospital 1
sd(theta.rank[1,])    # standard deviation of hospital 2
mean(theta.rank[1,] < n/2)  # probability hospital 1 falls in lower half of the hospital ranks 

# histograms, traces, autocorrelations for alpha and beta
par(mfrow=c(3,2))
# histograms 
hist(alpha.b)
hist(beta.b)
#traces
plot(alpha,pch=".")
plot(beta,pch=".")
# autocorrelations
acf(alpha)
acf(beta)
# mean ordered rank for all hospitals
mean.theta.rank = c()
for(i in 1:n) {
	mean.theta.rank[i] = mean(theta.rank[i,])
}
mean.rank = round(mean.theta.rank)
hospital.number = seq(1,n,1)
number.rank = cbind(hospital.number,mean.rank)
number.rank 



