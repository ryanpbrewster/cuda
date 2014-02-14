// pi1.cu
/*
 * A simple CUDA-enabled program that approximates \pi using monte-carlo
 * sampling. This version generates all the random numbers at the start,
 * then launches kernels to use them.
 */

#include <stdio.h>
#include <curand.h>

__global__ void pi(float* d_out, float* d_rands, int rands_per_kernel, int trials) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    int rand_idx = idx*rands_per_kernel;
    int count = 0;
    for(int i=0; i < trials; i++) {
        float x = d_rands[rand_idx + 2*i];
        float y = d_rands[rand_idx + 2*i + 1];
        if( x*x + y*y <= 1.0f ) {
            count++;
        }
    }
    d_out[idx] = float(count)/float(trials);
}

int main() {
    const int N_BLOCKS = 1024;
    const int N_THREADS = 512;
    const int N_KERNELS = N_BLOCKS * N_THREADS;
    const int N_TRIALS = 100;
    const int N_RANDS_PER_TRIAL = 2;
    const int N_RANDS = N_KERNELS * N_TRIALS * N_RANDS_PER_TRIAL;

    float* d_pis;
    float* d_rands;

    cudaMalloc(&d_pis, N_KERNELS * sizeof(float));
    cudaMalloc(&d_rands, N_RANDS * sizeof(float));

    curandGenerator_t prng;
    curandCreateGenerator(&prng, CURAND_RNG_PSEUDO_DEFAULT);
    curandSetPseudoRandomGeneratorSeed(prng, 0);
    curandGenerateUniform(prng, d_rands, N_RANDS);
    curandDestroyGenerator(prng);

    pi<<<N_BLOCKS, N_THREADS>>>(d_pis, d_rands, N_TRIALS*N_RANDS_PER_TRIAL, N_TRIALS);

    float* h_pis = (float*) malloc(N_KERNELS*sizeof(float));
    cudaMemcpy(h_pis, d_pis, N_KERNELS*sizeof(float), cudaMemcpyDeviceToHost);

    float avg = 0.0;
    for(int i=0; i < N_KERNELS; i++) {
        avg += h_pis[i];
    }
    avg /= N_KERNELS;

    printf("pi = %f\n", 4*avg);

    free(h_pis);
    cudaFree(d_pis);
    cudaFree(d_rands);
    return 0;
}
