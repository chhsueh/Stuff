########
library(MASS)

log_f = function(m=m,y=y,X=X,beta,beta.0,Sigma.0.inv){
	(sum(y*(X%*%beta)-m*(log(1+exp(X%*%beta))))-0.5*(t(beta-beta.0)%*%Sigma.0.inv%*%(beta-beta.0)))
}


## InI: starting point
bayes.logreg <- function(m,y,X,InI,beta.0,Sigma.0.inv=solve(1000*diag(11)),niter=10000,burnin=1000,
				 print.every=1000,retune=100,verbose=TRUE){
	
	## The we choose the initial proposal mean and var-covariance matrix by using the prior's.
	theta = InI 
	Sigma = diag(11)
	accept = 0
	accept_par = c()
	candidate_par = c()	
	par = c()

	for(kk in 1:(niter+burnin)){			
		if(kk <= burnin){

			## sample 1--4 ##
			theta_con = as.matrix(mvrnorm(1,theta,Sigma))
			rsample = theta_con[1:4]
			theta_t = as.matrix(c(rsample,theta[5:11]))
			u = log(runif(1))
			alpha = log_f(m=m,y=y,X=X,theta_t,beta.0,Sigma.0.inv) -
				   log_f(m=m,y=y,X=X,theta,beta.0,Sigma.0.inv)
			if (u < alpha){
				theta = theta_t 
				accept = accept+1
				accept_par = cbind(accept_par,theta)}
			par = cbind(par,theta)
			
			## sample 5--8 ##
			theta_con = as.matrix(mvrnorm(1,theta,Sigma))
			rsample = theta_con[5:8]
			theta_t = as.matrix(c(theta[1:4],rsample,theta[9:11]))
			u = log(runif(1))
			alpha = log_f(m=m,y=y,X=X,theta_t,beta.0,Sigma.0.inv) -
				   log_f(m=m,y=y,X=X,theta,beta.0,Sigma.0.inv)
			if (u < alpha){
				theta = theta_t 
				accept = accept+1
				accept_par = cbind(accept_par,theta)}
			par = cbind(par,theta)

			## sample 9--11 ##
			theta_con = as.matrix(mvrnorm(1,theta,Sigma))
			rsample = theta_con[9:11]
			theta_t = as.matrix(c(theta[1:8],rsample))
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
				r = accept/(3*retune)
				if(accept > 1){Sigma = (1-r)*Sigma+r*cov(t(accept_par))}
				if(accept < 3*retune/5){Sigma = 0.5*Sigma}
				print(list(paste("accept_ratio=",r,";","beta1=",theta[1],";","beta2=",theta[2],sep=""),Sigma=Sigma))
				accept_par=c()
				accept = 0}

		}else {

			## sample 1--4 ##
			theta_con = as.matrix(mvrnorm(1,theta,Sigma))
			rsample = theta_con[1:4]
			theta_t = as.matrix(c(rsample,theta[5:11]))
			u = log(runif(1))
			alpha = log_f(m=m,y=y,X=X,theta_t,beta.0,Sigma.0.inv) -
				   log_f(m=m,y=y,X=X,theta,beta.0,Sigma.0.inv)
			if (u < alpha){theta = theta_t; accept = accept+1}
			candidate_par = cbind(candidate_par,theta)
			par = cbind(par,theta)
			
			## sample 5--8 ##
			theta_con = as.matrix(mvrnorm(1,theta,Sigma))
			rsample = theta_con[5:8]
			theta_t = as.matrix(c(theta[1:4],rsample,theta[9:11]))
			u = log(runif(1))
			alpha = log_f(m=m,y=y,X=X,theta_t,beta.0,Sigma.0.inv) -
				   log_f(m=m,y=y,X=X,theta,beta.0,Sigma.0.inv)
			if (u < alpha){theta = theta_t; accept = accept+1}
			candidate_par = cbind(candidate_par,theta)
			par = cbind(par,theta)

			## sample 9--11 ##
			theta_con = as.matrix(mvrnorm(1,theta,Sigma))
			rsample = theta_con[9:11]
			theta_t = as.matrix(c(theta[1:8],rsample))
			u = log(runif(1))
			alpha = log_f(m=m,y=y,X=X,theta_t,beta.0,Sigma.0.inv) -
				   log_f(m=m,y=y,X=X,theta,beta.0,Sigma.0.inv)
			if (u < alpha){theta = theta_t; accept = accept+1}
			candidate_par = cbind(candidate_par,theta)
			par = cbind(par,theta)


			## Print an update to the user after every period of this many iterations.
			if (kk%%print.every == 0){
				print(paste("accept_rate=",accept/(3*print.every),";","beta1=",theta[1],";","beta2=",theta[2],sep=""));accept=0}

		} # end of if loop

	} # end of for loop

	return=list(Can=candidate_par)
	

} # end of function 

###################################################################
#######################  Fit the model  ###########################
###################################################################

b = matrix(0,11,1)
Data = read.table("breast_cancer.txt",header=TRUE)
M = which(Data$diagnosis=="M")
Data$y = matrix(0,dim(Data)[1],1)
Data$y[M] = 1

m = as.matrix(rep(1,dim(Data)[1]))
y = as.matrix(Data$y)
X = cbind(rep(1,dim(Data)[1]),as.matrix(Data[,1:10]))

### Standardlize the covariates of X ####
for (i in 2:ncol(X)){
	m.X = mean(X[,i]); sd.X = sd(X[,i])
	X[,i] = (X[,i]-m.X)/sd.X
}


beta.0 = matrix(0,11,1)
Sigma.0.inv = solve(1000*diag(11))

Get.initial =  bayes.logreg(m,y,X,beta.0,beta.0,Sigma.0.inv,niter=20000,burnin=5000,print.every=5000,retune=100,verbose=TRUE)
Ini = as.matrix(Get.initial$Can[,length(Get.initial$Can[1,])])
Bayes = bayes.logreg(m,y,X,Ini,beta.0,Sigma.0.inv,niter=20000,burnin=3000,print.every=5000,retune=100,verbose=TRUE)
Beta = Bayes$Can

filename <- paste0("prob3.csv")
write.table(data.frame(Beta),file=filename,sep=",",quote=FALSE,row.names=FALSE,col.names=FALSE)


### Check the algorithm (Convergence) ###

for (kk in 1:100){

	## step 1 ##
	theta = mvrnorm(1,as.matrix(rep(0,11)),1000*diag(11))

	## step 2 ##
	logit.inv = function(x){exp(x)/(1+exp(x))}
	Y = apply(X, 1, function(x) rbinom(1, 1, logit.inv(x%*%as.matrix(theta))) )
	Y[is.na(Y)] = 0

	## step 3 ##
	Initial.pt = as.matrix(theta)
	Bayes = bayes.logreg(m,as.matrix(Y),X,beta.0,beta.0,Sigma.0.inv,niter=20000,burnin=3000,print.every=5000,retune=100,verbose=TRUE)

	CI = c()
	for (kk in 1:11){
	q = quantile(Bayes$Can[kk,],c(0.025,0.975))
	CI = rbind(CI,q)
	}	

	## step 4 ##
	for (i in 1:11){
	if (findInterval(theta[i],CI[i,])==1){b[i]=b[i]+1}}
		
}

niter=20000;burnin=3000
check = c(b,niter,burnin)
names(check)= c(as.character(1:11),"niter","burnin")

filename <- paste0("Check_MHGIB.csv")
write.table(data.frame(check),file=filename,sep=",",quote=FALSE,row.names=TRUE,col.names=FALSE)


#######################################################################################################
Prob3 = read.csv("prob3.csv",header=FALSE)

## c ##
CI = c()
for (kk in 1:11){
	q = quantile(Prob3[kk,],c(0.025,0.975))
	CI = rbind(CI,q)
}
CI = as.data.frame(CI)
row.names(CI)=c("intercept","area" ,"compactness", "concavepts", "concavity", "fracdim", "perimeter", "radius", "smoothness", "symmetry", "texture")

## acf function ##
Auto_cor = c()
for (i in 1:11){
	ac = acf(t(Prob3[i,]),lag=1,plot=FALSE)
	Auto_cor = rbind(Auto_cor,ac$acf[2])}

## predictive checking ##
## step 1 ##theta = as.matrix(Prob3[,8000:10000])

## step 2 ##
s = c()
for (i in 1:2001){
	logit.inv = function(x){exp(x)/(1+exp(x))}
	Y = apply(X, 1, function(x) rbinom(1, 1, logit.inv(x%*%as.matrix(theta[,i]))) )
	s = c(s,mean(Y))
}

	
## plot convergence ##
par(mfrow=c(4,3))
for (i in 1:11){
	plot(1:10000,Prob3[i,],type="l",ylab = "beta",main = paste("beta",i,sep="."))
}







