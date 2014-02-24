// curand2.cu
/*
 * A simple CUDA-enabled program that generates random numbers on-the-fly
 * within each kernel.
 */

#include <stdio.h>
#include <curand.h>
#include <curand_kernel.h>

__global__ void rnd(curandState_t* states) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    curand_init(0, idx, 0, &states[idx]);

    printf("Thread (%d,%d) --> %f\n", blockIdx.x, threadIdx.x, curand_uniform(&states[idx]));
}

int main() {
    const int N_BLOCKS = 1;
    const int N_THREADS = 16;
    const int N_KERNELS = N_BLOCKS * N_THREADS;

    curandState_t* states;
    cudaMalloc(&states, N_KERNELS * sizeof(curandState_t));

    rnd<<<N_BLOCKS, N_THREADS>>>(states);

    cudaFree(states);
    return 0;
}
