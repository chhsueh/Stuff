library(truncnorm)
library(MASS)
library(compiler)
library(RCUDA)

m = loadModule("rtruncnorm.ptx")
my_kernel <- m$rtruncnorm_kernel

"compute_grid" <- function(N,sqrt_threads_per_block=16L,grid_nd=1)
{
    # if...
    # N = 1,000,000
    # => 1954 blocks of 512 threads will suffice
    # => (62 x 32) grid, (512 x 1 x 1) blocks
    # Fix block dims:
    block_dims <- c(as.integer(sqrt_threads_per_block), as.integer(sqrt_threads_per_block), 1L)
    threads_per_block <- prod(block_dims)
    if (grid_nd==1){
      grid_d1 <- as.integer(max(1L,ceiling(N/threads_per_block)))
      grid_d2 <- 1L
    } else {
      grid_d1 <- as.integer(max(1L, floor(sqrt(N/threads_per_block))))
      grid_d2 <- as.integer(ceiling(N/(grid_d1*threads_per_block)))
    }
    grid_dims <- c(grid_d1, grid_d2, 1L)
    return(list("grid_dims"=grid_dims,"block_dims"=block_dims))
}

z_sample  = function(gpu = FALSE, y, X, beta, rngnum){
	Z = matrix(0,length(y),1)
	if (gpu){
		#compute the grid
		n = as.integer(length(y))
		grid = compute_grid(n,sqrt_threads_per_block=16L,grid_nd=1)
		grid_dims <- as.integer(grid$grid_dims)
		block_dims <- as.integer(grid$block_dims)	

		x = matrix(0,n,1)
		mu = X%*%beta
		sigma = matrix(1,n,1)
		w1 = which(y==1)
		w2 = which(y==0)
		lo=matrix(0,n,1)
		lo[w2] = as.integer(-10000) 
		hi=matrix(0,n,1)
		hi[w1] = as.integer(10000)
		maxtries = 5000L
	
		Z= .cuda(my_kernel,"x"=x,n,mu,sigma,lo,hi,maxtries,rngnum,gridDim=grid_dims,blockDim=block_dims,outputs="x")
		Z = as.matrix(Z)

	} else{
		for(i in 1:length(y)){
		if(y[i] == 0){
			Z[i] = rtruncnorm(1, a=-Inf, b=0, mean = X[i,]%*%beta, sd = 1) 
		} else{
			Z[i] = rtruncnorm(1, a=0, b=Inf, mean = X[i,]%*%beta, sd = 1)
		}	
		}
	}
	return(Z)	
}
z_sample = cmpfun(z_sample)

probit_mcmc = function(
      y,           # vector of length n 
      X,           # (n x p) design matrix
      beta_0,      # (p x 1) prior mean
      Sigma_0_inv, # (p x p) prior precision 
      niter = 2000,       # number of post burnin iterations
      burnin = 500,      # number of burnin iterations
      gpu = FALSE # specify using gpu or cpu
      ){
	
	##initial value
	X = as.matrix(X)
	p = ncol(X)
	B = solve(Sigma_0_inv+t(X)%*%X)
	Z = matrix(0,length(y),1)
	beta = beta_0
	Beta = c()

	for(kk in 1:(niter+burnin)){
		seed = set.seed(kk)
		rngnum = as.integer(kk)			
		if(kk <= burnin){

			Z = z_sample(gpu, y, X, beta, rngnum)
			beta = as.matrix(mvrnorm(1,B%*%(Sigma_0_inv%*%beta_0+t(X)%*%Z),B))

		} else{

			Z = z_sample(gpu, y, X, beta, rngnum)
			beta = as.matrix(mvrnorm(1,B%*%(Sigma_0_inv%*%beta_0+t(X)%*%Z),B))
			Beta = cbind(Beta,beta)

		} # end of if
	} # end of for loop
	return(Beta)

} # end of fuction
probit_mcmc = cmpfun(probit_mcmc)

### Try mini-data ###
Data = read.table("mini_data.txt",header=TRUE)
X = Data[,2:ncol(Data)]
y = as.matrix(Data$y)
p = ncol(X)
beta_0 = matrix(0,p,1)
Sigma_0_inv = matrix(0,p,p)
niter = 2000
burnin =500

Systime = matrix(0,2,5)

List <- list.files(pattern = c("data_0"), all.files =FALSE) 
for (i in 1:3){
	
	Data = read.table(List[i],header=TRUE)
	X = Data[,2:ncol(Data)]
	y = as.matrix(Data$y)
	p = ncol(X)
	beta_0 = matrix(0,p,1)
	Sigma_0_inv = matrix(0,p,p)

	T1=system.time(probit_mcmc(y, X, beta_0, Sigma_0_inv, niter=2000, burnin=500, gpu = FALSE))
	T2=system.time(probit_mcmc(y, X, beta_0, Sigma_0_inv, niter, burnin, gpu = TRUE))
	Systime[1,i] = T1[3]
	Systime[2,i] = T2[3]

}

save(Systime,file="MCMC_time.RData")

i = 4
	
	Data = read.table(List[i],header=TRUE)
	X = Data[,2:ncol(Data)]
	y = as.matrix(Data$y)
	p = ncol(X)
	beta_0 = matrix(0,p,1)
	Sigma_0_inv = matrix(0,p,p)

	T1=system.time(probit_mcmc(y, X, beta_0, Sigma_0_inv, niter=2000, burnin=500, gpu = FALSE))
	T2=system.time(probit_mcmc(y, X, beta_0, Sigma_0_inv, niter, burnin, gpu = TRUE))
	Systime[1,i] = T1[3]
	Systime[2,i] = T2[3]

save(Systime,file="MCMC_time2.RData")




#####
#z_sample(gpu = TRUE, y, X, beta, rngnum)
