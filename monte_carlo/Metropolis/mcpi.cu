// metropolis1.cu
/*
 * A simple CUDA-enabled program that approximates Pi by evaluating
 *     Integrate[ Sqrt[1-x^2], {x,-1,1} ]
 * using Metropolis Monte Carlo, with the weight function A*(1-x^2)
 * where A is a normalization factor
 */

#include <iostream>
#include <curand.h>
#include <curand_kernel.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdbool.h>

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}


__device__ double weightFunction(double x) {
    return 0.75*(1-x*x);
}
__device__ double f(double x) {
    return sqrt(1-x*x);
}





__global__ void initThreads(double* d_out, curandState_t* states, double* xs, double* radii, int seed) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    curand_init(seed+idx, 0, 0, &states[idx]);
    d_out[idx] = 0.0;
    xs[idx] = 0.0;
    radii[idx] = 1.0;
}


__global__ void mcpi(double* d_out, curandState_t* states, double* xs, double* radii, int N_TRIALS) {
    const int ITERS_PER_RADIUS_ADJ = 100;

    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    curandState_t tmp_state = states[idx];
    double x = xs[idx];
    double r = radii[idx];

    int accept = 0, reject = 0;
    double estimate = 0.0;
    double cur_weight = weightFunction(x);
    for(int i=1; i <= N_TRIALS; i++) {
        double xp = x + r*(2*curand_uniform_double(&tmp_state)-1);
        double new_weight = weightFunction(xp);
        double xi = curand_uniform_double(&tmp_state);
        if( xi < new_weight/cur_weight ) {
            x = xp;
            cur_weight = new_weight;
            accept++;
        } else {
            reject++;
        }
        estimate += f(x) / cur_weight;

        if( i%ITERS_PER_RADIUS_ADJ == 0 ) {
            double accept_ratio = double(accept) / double(accept + reject);
            double adj = min(max(2.0*accept_ratio, 0.9), 1.1);
            r *= adj;
            accept = 0; reject = 0;
        }
    }

    states[idx] = tmp_state;
    xs[idx] = x;
    radii[idx] = r;
    d_out[idx] += estimate / N_TRIALS;
}

double pi(int const GRID_SIZE, int const BLOCK_SIZE, int const N_THREADS, int const N_RUNS, int seed) {
    int N_KERNELS = GRID_SIZE * BLOCK_SIZE;

    double* h_pis = (double*) malloc(N_KERNELS*sizeof(double));
    double* d_pis;
    gpuErrchk(cudaMalloc(&d_pis, N_KERNELS * sizeof(double)));

    curandState_t* states;
    gpuErrchk(cudaMalloc(&states, N_KERNELS * sizeof(curandState_t)));

    double* xs;
    gpuErrchk(cudaMalloc(&xs, N_KERNELS * sizeof(double)));

    double* radii;
    gpuErrchk(cudaMalloc(&radii, N_KERNELS * sizeof(double)));

    initThreads<<<GRID_SIZE, BLOCK_SIZE>>>(d_pis, states, xs, radii, seed);
    for(int irun=1; irun <= N_RUNS; irun++) {
        pi<<<GRID_SIZE, BLOCK_SIZE>>>(d_pis, states, xs, radii, N_TRIALS);
    }
    gpuErrchk(cudaMemcpy(h_pis, d_pis, N_KERNELS*sizeof(double), cudaMemcpyDeviceToHost));

    double avg = 0.0;
    for(int i=0; i < N_KERNELS; i++) {
        avg += h_pis[i] / N_RUNS;
    }
    avg /= N_KERNELS;

    free(h_pis);
    cudaFree(d_pis);
    cudaFree(states);
    cudaFree(xs);
    return avg;
}
