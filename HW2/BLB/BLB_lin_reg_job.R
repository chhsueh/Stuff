
mini <- FALSE

#============================== Setup for running on Gauss... ==============================#

args <- commandArgs(TRUE)

cat("Command-line arguments:\n")
print(args)

####
# sim_start ==> Lowest possible dataset number
###

###################
sim_start <- 1000
###################

if (length(args)==0){
  sim_num <- sim_start + 1
  set.seed(121231)
} else {
  # SLURM can use either 0- or 1-indexing...
  # Lets use 1-indexing here...
  sim_num <- sim_start + as.numeric(args[1])
  sim_seed <- (762*(sim_num-1) + 121231)
}

cat(paste("\nAnalyzing dataset number ",sim_num,"...\n\n",sep=""))

# Find r and s indices:
num = sim_num-1000
if (num%%50==0){
	r_index = 50
	s_index = floor(num/50)
}else{
	r_index = num%%50
	s_index = floor(num/50)+1

}

#============================== Run the simulation study ==============================#

# Load packages:
library(BH)
library(bigmemory.sri)
library(bigmemory)
library(biganalytics)

# I/O specifications:
datapath <- "/home/pdbaines/data"
outpath <- "output/"

# mini or full?
if (mini){
	rootfilename <- "blb_lin_reg_mini"
} else {
	rootfilename <- "blb_lin_reg_data"
}

# Filenames:
filename = paste(rootfilename,".desc",sep="")

# Set up I/O stuff:
file = file.path(datapath,filename)

# Attach big.matrix :
Data = attach.big.matrix(file)

# Remaining BLB specs:
ga = 0.7
n = dim(Data)[1]
b = floor(n^ga)

# Extract the subset:
set.seed(s_index)
S = sample(1:n, b, replace=FALSE)

# Reset simulation seed:
rm(.Random.seed)

# Bootstrap dataset:
boot.sample = rmultinom(1,n,prob=rep(1/b,b))

# Fit lm:
fit = lm(Data[S,dim(Data)[2]]~Data[S,1:(dim(Data)[2]-1)]-1,weights = boot.sample)
coeff = fit$coefficient

# Output file:
outfile = paste0("output/","coef_",sprintf("%02d",s_index),"_",sprintf("%02d",r_index),".txt")

# Save estimates to file:
write.table(coeff,file=outfile,row.names=FALSE)



