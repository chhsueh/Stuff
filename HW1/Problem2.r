############################        PROBLEM2         ######################################
###########################################################################################

Data = read.csv("blr_data_1001.csv",header=TRUE)
m = as.matrix(Data$n)
y = as.matrix(Data$y)
X = as.matrix(Data[,3:4])
library(MASS)

## 2-d ##
log_f = function(m=m,y=y,X=X,beta,beta.0=as.matrix(c(0,0)),Sigma.0.inv=diag(2)){
	(sum(y*(X%*%beta)-m*(log(1+exp(X%*%beta))))-0.5*(t(beta-beta.0)%*%Sigma.0.inv%*%(beta-beta.0)))
}
#f = function(m=m,y=y,X=X,beta,beta.0=as.matrix(c(0,0)),Sigma.0.inv=diag(2)){
#	exp((sum(y*(X%*%beta)-m*(log(1+exp(X%*%beta))))-0.5*(t(beta-beta.0)%*%Sigma.0.inv%*%(beta-beta.0))))
#}


bayes.logreg <- function(m,y,X,beta.0=as.matrix(c(0,0)),Sigma.0.inv=diag(2),niter=10000,burnin=1000,
				 print.every=1000,retune=100,verbose=TRUE){
	
	## The we choose the initial proposal mean and var-covariance matrix by using the prior's.
	theta = beta.0
	Sigma = solve(Sigma.0.inv)
	accept = 0
	accept_par = c()
	candidate_par = c()	
	par = c()

	for(kk in 1:(niter+burnin)){			
		if(kk <= burnin){

			theta_t = as.matrix(mvrnorm(1,theta,Sigma))
			u = log(runif(1))
			alpha = log_f(m=m,y=y,X=X,theta_t,beta.0,Sigma.0.inv) -
				   log_f(m=m,y=y,X=X,theta,beta.0,Sigma.0.inv)
			if (u < alpha){
				theta = theta_t 
				accept = accept+1
				accept_par = cbind(accept_par,theta)}
			par = cbind(par,theta)

			## retune the variance of the proposal distribution every retune iterations.
			if (kk%%retune == 0){
				r = accept/retune
				if(accept > 1){Sigma = (1-r)*Sigma+r*cov(t(accept_par))}
				if(accept < retune/3){Sigma = 0.5*Sigma}
				print(list(paste("accept_ratio=",r,";","beta1=",theta[1],";","beta2=",theta[2],sep=""),Sigma=Sigma))
				accept_par=c()
				accept = 0}

		}else {

			theta_t = as.matrix(mvrnorm(1,theta,Sigma))
			u = log(runif(1))
			alpha = log_f(m=m,y=y,X=X,theta_t,beta.0,Sigma.0.inv) -
				   log_f(m=m,y=y,X=X,theta,beta.0,Sigma.0.inv)
			if (u < alpha){theta = theta_t; accept = accept+1}
			candidate_par = cbind(candidate_par,theta)
			par = cbind(par,theta)

			## Print an update to the user after every period of this many iterations.
			if (kk%%print.every == 0){
				print(paste("accept_rate=",accept/print.every,";","beta1=",theta[1],";","beta2=",theta[2],sep=""));accept=0}

		} # end of if loop

	} # end of for loop

	Quantile = rbind(q.beta1 = quantile(candidate_par[1,],seq(0.01,0.99,0.01)),
			     quantile(candidate_par[2,],seq(0.01,0.99,0.01)))

	
	C.I.beta1 = quantile(candidate_par[1,],c(0.025,0.975))
	C.I.beta2 = quantile(candidate_par[2,],c(0.025,0.975))
	return=list(C.I.beta1=C.I.beta1, C.I.beta2=C.I.beta2, Quantile=Quantile)
	

} # end of function 

bayes.logreg(m,y,X,beta.0=as.matrix(c(0,0)),Sigma.0.inv=diag(2),niter=10000,burnin=1000,print.every=1000,retune=100,verbose=TRUE)     


##############################################################################
##      1,2,...,99% percentiles of the marginal posterior distributions     ##
##############################################################################

q.beta1 = quantile(candidate_par[1,],seq(0.01,0.99,0.01))
q.beta2 = quantile(candidate_par[2,],seq(0.01,0.99,0.01))


###############################################################################
##             Simulated study to check everything is working                ##
###############################################################################
Parameter = c()
Cre.int = c()
b1 = 0
b2 = 0
both = 0

for (kk in 1:200){

	## step 1 ##
	theta = mvrnorm(1,as.matrix(c(0,0)),diag(2))
	Parameter = rbind(Parameter,theta)

	## step 2 ##
	logit.inv = function(x){exp(x)/(1+exp(x))}
	Y = apply(Data, 1, function(x) sum(rbinom(x[2], 1, logit.inv(x[3:4]%*%as.matrix(theta)))) )

	## step 3 ##
	Bayes = bayes.logreg(m,as.matrix(Y),X,as.matrix(theta),Sigma.0.inv=diag(2),niter=10000,burnin=1000,print.every=1000,retune=100,verbose=TRUE)

	## step 4 ##
	if (findInterval(theta[1],Bayes$C.I.beta1)==1){b1=b1+1}
	if (findInterval(theta[2],Bayes$C.I.beta2)==1){b2=b2+1}
	if (findInterval(theta[1],Bayes$C.I.beta1)==1 & findInterval(theta[2],Bayes$C.I.beta2)==1){both=both+1}
	
}






