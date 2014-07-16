// pi2.cu
/*
 * A simple CUDA-enabled program that approximates \pi using monte-carlo
 * sampling. This version generates random numbers on-the-fly within each
 * kernel.
 */

#include <iostream>
#include <curand.h>
#include <curand_kernel.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdbool.h>
using namespace std;

__global__ void initThreads(float* d_out, curandState_t* states) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    curand_init(idx, idx, 0, &states[idx]);
    d_out[idx] = 0.0;
}

__global__ void pi(float* d_out, curandState_t* states, int N_TRIALS) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    int count = 0;
    for(int i=1; i <= N_TRIALS; i++) {
        float x = curand_uniform(&states[idx]);
        float y = curand_uniform(&states[idx]);
        if( x*x + y*y <= 1.0f ) {
            count++;
        }
    }

    d_out[idx] += float(count)/float(N_TRIALS);
}

int main(int argc, char** argv) {
    int GRID_SIZE  = 256;
    int BLOCK_SIZE = 256;
    int N_TRIALS   = 1000;
    int N_RUNS     = 10;

    char x;
    opterr = 0;
    while((x = getopt(argc, argv, "g:b:t:r:")) != -1) {
        switch(x) {
            case 'g': GRID_SIZE  = atoi(optarg); break;
            case 'b': BLOCK_SIZE = atoi(optarg); break;
            case 't': N_TRIALS   = atoi(optarg); break;
            case 'r': N_RUNS     = atoi(optarg); break;
            case '?':
                if (optopt == 'g' || optopt == 'b' || optopt == 't' || optopt == 'r') {
                    fprintf (stderr, "Option -%c requires an argument.\n", optopt);
                } else {
                    fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                }
                abort();
            default: abort();
        }
    }

    int N_KERNELS = GRID_SIZE * BLOCK_SIZE;

    float* h_pis = (float*) malloc(N_KERNELS*sizeof(float));
    float* d_pis;
    cudaMalloc(&d_pis, N_KERNELS * sizeof(float));

    curandState_t* states;
    cudaMalloc(&states, N_KERNELS * sizeof(curandState_t));

    time_t start = clock();

    initThreads<<<GRID_SIZE, BLOCK_SIZE>>>(d_pis, states);
    for(int irun=1; irun <= N_RUNS; irun++) {
        pi<<<GRID_SIZE, BLOCK_SIZE>>>(d_pis, states, N_TRIALS);
    }
    cudaMemcpy(h_pis, d_pis, N_KERNELS*sizeof(float), cudaMemcpyDeviceToHost);
    float avg = 0.0;
    for(int i=0; i < N_KERNELS; i++) {
        avg += h_pis[i] / N_RUNS;
    }
    avg /= N_KERNELS;

    time_t end = clock();

    int64_t iters = int64_t(N_KERNELS)*int64_t(N_TRIALS)*int64_t(N_RUNS);
    int elapsed = 1000*(end-start)/CLOCKS_PER_SEC;

    cout << float(iters)/float(elapsed) << " iters/ms\n";

    free(h_pis);
    cudaFree(d_pis);
    cudaFree(states);
    return 0;
}
