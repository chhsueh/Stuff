
##
#
# Logistic regression
# 
# Y_{i} | \beta \sim \textrm{Bin}\left(n_{i},e^{x_{i}^{T}\beta}/(1+e^{x_{i}^{T}\beta})\right)
# \beta \sim N\left(\beta_{0},\Sigma_{0}\right)
#
##

library(mvtnorm)
library(coda)
library(MASS)

########################################################################################
########################################################################################
## Handle batch job arguments:

# 1-indexed version is used now.
args <- commandArgs(TRUE)

cat(paste0("Command-line arguments:\n"))
print(args)

####
# sim_start ==> Lowest simulation number to be analyzed by this particular batch job
###

#######################
sim_start <- 1000
length.datasets <- 200
#######################

if (length(args)==0){
  sinkit <- FALSE
  sim_num <- sim_start + 1
  set.seed(1330931)
} else {
  # Sink output to file?
  sinkit <- TRUE
  # Decide on the job number, usually start at 1000:
  sim_num <- sim_start + as.numeric(args[1])
  # Set a different random seed for every job number!!!
  set.seed(762*sim_num + 1330931)
}

# Simulation datasets numbered 1001-1200

########################################################################################
########################################################################################

log_f <- function(m=m,y=y,X=X,beta,beta.0=as.matrix(c(0,0)),Sigma.0.inv=diag(2)){
	(sum(y*(X%*%beta)-m*(log(1+exp(X%*%beta))))-0.5*(t(beta-beta.0)%*%Sigma.0.inv%*%(beta-beta.0)))}


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
			alpha = log_f(m=m,y=y,X=X,theta_t,beta.0,Sigma.0.inv=diag(2)) -
				   log_f(m=m,y=y,X=X,theta,beta.0,Sigma.0.inv=diag(2))
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
			alpha = log_f(m=m,y=y,X=X,theta_t,beta.0,Sigma.0.inv=diag(2)) -
				   log_f(m=m,y=y,X=X,theta,beta.0,Sigma.0.inv=diag(2))
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

#################################################
# Set up the specifications:
beta.0 <- matrix(c(0,0))
Sigma.0.inv <- diag(rep(1.0,2))
niter <- 10000
# etc... (more needed here)
#################################################

DIR = "data"
file = paste("blr_data_",sim_num,".csv",sep="")
f.p <- file.path(DIR, file)

# Read data corresponding to appropriate sim_num:
Data = read.csv(f.p, header=TRUE)

# Extract X and y:
m = as.matrix(Data$n)
y = as.matrix(Data$y)
X = as.matrix(Data[,3:4])

# Fit the Bayesian model:
Bayes = bayes.logreg(m,y,X,beta.0,Sigma.0.inv=diag(2),niter=10000,burnin=1000,print.every=1000,retune=100,verbose=TRUE)

# Extract posterior quantiles...
Quantile = t(Bayes$Quantile)

# Write results to a (99 x p) csv file...
dir <- "results/"
filename <- paste0(dir,"blr_res_",sim_num,".csv")
write.table(data.frame(Quantile),file=filename,sep=",",quote=FALSE,row.names=FALSE,col.names=FALSE)

# Go celebrate.
 
cat("done. :)\n")







