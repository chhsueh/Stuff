#include <stdio.h>
#include <stdlib.h>

#include <cuda.h>
#include <curand_kernel.h>
#include <math_constants.h>

extern "C"
{

__global__ void 
rtruncnorm_kernel(
    float *x, 
    int n, 
    float *mu, 
    float *sigma, 
    float *lo, 
    float *hi,
    int maxtries,
    int rngnum) //number of the random seed
{
    int accepted;
    float sample;
    int numtries;
    float m;
    float alpha;

    //variables for rejection sampling
    float rexp; 
    float z; 
    float phi;
    float u;

    // Usual block/thread indexing...
    int myblock = blockIdx.x + blockIdx.y * gridDim.x;
    int blocksize = blockDim.x * blockDim.y * blockDim.z;
    int subthread = threadIdx.z*(blockDim.x * blockDim.y) + threadIdx.y*blockDim.x + threadIdx.x;
    
    // map these into a single index : idx
    int idx = myblock * blocksize + subthread;

    // check whether the idx is <n
    if (idx < n)
    {

        // Setup the RNG:
        curandState rng_state;
        curand_init(9131+idx*17,rngnum,0,&rng_state);

        // Sample:
        accepted = 0;
        numtries = 0;

        while (accepted == 0 && numtries < maxtries)
        {
            sample = mu[idx]+sigma[idx]*curand_normal(&rng_state);
            numtries = numtries+1;

            if (sample>lo[idx] && sample<hi[idx])
            {

                accepted = 1;
                x[idx] = sample;
                //printf("rnumber = %f\n", sample);

            } //end of if(small) loop

        } // end of while loop

	while (accepted == 0) //if accepted = 0 run rejection sampling.
	{ 
	    //code for rejection sampling.
	    if(abs(lo[idx]-mu[idx]) < abs(hi[idx]-mu[idx])){ //right tail

		m = abs((lo[idx]-mu[idx])/sigma[idx]);
		alpha = (m+sqrt(pow(m,2)+4))/2;
	    	rexp = -log(curand_uniform(&rng_state))/alpha;
	    	z = m + rexp;
		if (m<alpha){
			phi = exp(-pow(alpha-z,2)/2);
	    	}
	    	else{
			phi = exp(pow(m-alpha,2)/2-pow(alpha-z,2)/2);
	    	} //decide phi

		u = curand_uniform(&rng_state);
		if (u<phi){
			accepted = 1;
			x[idx] = mu[idx]+sigma[idx]*z;
	    	} 
		
	    } else{ //left tail

		m = abs((mu[idx]-hi[idx])/sigma[idx]);
		alpha = (m+sqrt(pow(m,2)+4))/2;
	    	rexp = -log(curand_uniform(&rng_state))/alpha;
	    	z = m + rexp;
		if (m<alpha){
		phi = exp(-pow(alpha-z,2)/2);
	    	}
	    	else{
			phi = exp(pow(m-alpha,2)/2-pow(alpha-z,2)/2);
	    	} //decide phi

		u = curand_uniform(&rng_state);
		if (u<phi){
			accepted = 1;
			x[idx] = mu[idx]-sigma[idx]*z;
	    	} 

	    }
	
		
	} // end of rejection sampling.

    } // end of if loop

    return;
} // end of function

} // END extern "C"


//#### More variables: ########################
                  //int mu_len, 
		  //int sigma_len,
                  //int lo_len, 
		  //int hi_len,
                  //int maxtries
//#############################################
