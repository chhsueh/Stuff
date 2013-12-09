
##### Function for computing default grid/block sizes #####


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


###### variables in truncnorm_kernel ################
    float *x, 
    int n, 
    float *mu, 
    float *sigma, 
    float *lo, 
    float *hi,
    int maxtries,
    int rngnum, //number of the random seed
#####################################################

## b ##

m = loadModule("rtruncnorm.ptx")
my_kernel <- m$rtruncnorm_kernel


n=as.integer(10000)
grid = compute_grid(n,sqrt_threads_per_block=16L,grid_nd=1)
grid_dims <- as.integer(grid$grid_dims)
block_dims <- as.integer(grid$block_dims)

x = matrix(0,n,1)
mu = matrix(2,n,1)
sigma = matrix(1,n,1)
lo=matrix(0,n,1)
hi=matrix(1.5,n,1)
maxtries = 5000L
rngnum = 27L
sample=.cuda(my_kernel,"x"=x,n,mu,sigma,lo,hi,maxtries,rngnum,gridDim=grid_dims,blockDim=block_dims,outputs="x")


## c ## comparing system time use system.time ##

library(truncnorm)
T = matrix(0,2,8)

for (i in 1:8)
{
	n=as.integer(10^i)
	grid = compute_grid(n,sqrt_threads_per_block=16L,grid_nd=1)

	x = matrix(0,n,1)
	mu = matrix(2,n,1)
	sigma = matrix(1,n,1)
	lo=matrix(0,n,1)
	hi=matrix(1.5,n,1)
	maxtries = 50L
	m = matrix(0,n,1)
	rngnum = 27L
	timeCPU = system.time(rtruncnorm(n, a=0, b=1.5, mean = 0, sd = 1))
	timeGPU = system.time({.cuda(my_kernel,"x"=x,n,mu,sigma,lo,hi,maxtries,rngnum,gridDim=grid_dims,blockDim=block_dims,outputs="x")})
	T[1,i] = timeCPU[3]
	T[2,i] = timeGPU[3]
}

T = as.data.frame(T)
rownames(T)=c("CPU","GPU")
names(T)=c(1:8)


#### g
n = as.integer(10000)
grid = compute_grid(n,sqrt_threads_per_block=16L,grid_nd=1)
grid_dims <- as.integer(grid$grid_dims)
block_dims <- as.integer(grid$block_dims)
x = matrix(0,n,1)
mu = matrix(0,n,1)
sigma = matrix(1,n,1)
lo=matrix(-10000,n,1)
hi=matrix(-10,n,1)
maxtries = 2000L
rngnum = 27L

sample_tail=.cuda(my_kernel,"x"=x,n,mu,sigma,lo,hi,maxtries,rngnum,gridDim=grid_dims,blockDim=block_dims,outputs="x")
