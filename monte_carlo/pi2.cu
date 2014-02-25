// pi2.cu
/*
 * A simple CUDA-enabled program that approximates \pi using monte-carlo
 * sampling. This version generates random numbers on-the-fly within each
 * kernel.
 */

#include <stdio.h>
#include <curand.h>
#include <curand_kernel.h>

__global__ void pi(double* d_out, curandState_t* states, int trials) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    curand_init(blockIdx.x, threadIdx.x, 0, &states[idx]);

    int count = 0;
    for(int i=0; i < trials; i++) {
        double x = curand_uniform_double(&states[idx]);
        double y = curand_uniform_double(&states[idx]);
        if( x*x + y*y <= 1.0f ) {
            count++;
        }
    }
    d_out[idx] = double(count)/double(trials);
}

int main() {
    const int N_BLOCKS = 1024;
    const int N_THREADS = 256;
    const int N_KERNELS = N_BLOCKS * N_THREADS;
    const int N_TRIALS = 10000;

    double* d_pis;
    cudaMalloc(&d_pis, N_KERNELS * sizeof(double));

    curandState_t* states;
    cudaMalloc(&states, N_KERNELS * sizeof(curandState_t));

    pi<<<N_BLOCKS, N_THREADS>>>(d_pis, states, N_TRIALS);

    double* h_pis = (double*) malloc(N_KERNELS*sizeof(double));
    cudaMemcpy(h_pis, d_pis, N_KERNELS*sizeof(double), cudaMemcpyDeviceToHost);

    double avg = 0.0;
    for(int i=0; i < N_KERNELS; i++) {
        avg += h_pis[i];
    }
    avg /= N_KERNELS;

    printf("pi = %f\n", 4*avg);

    free(h_pis);
    cudaFree(d_pis);
    cudaFree(states);
    return 0;
}
