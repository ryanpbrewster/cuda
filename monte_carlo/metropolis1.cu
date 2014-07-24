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
using namespace std;

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





__global__ void initThreads(double* d_out, curandState_t* states, double* xs, double* radii) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    curand_init(idx, idx, 0, &states[idx]);
    d_out[idx] = 0.0;
    xs[idx] = 0.0;
    radii[idx] = 1.0;
}


__global__ void pi(double* d_out, curandState_t* states, double* xs, double* radii, int N_TRIALS) {
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

    printf("G=% 4d, B=% 4d, T=% 8d, R=% 5d\n", GRID_SIZE, BLOCK_SIZE, N_TRIALS, N_RUNS);

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

    time_t start = clock();

    initThreads<<<GRID_SIZE, BLOCK_SIZE>>>(d_pis, states, xs, radii);
    for(int irun=1; irun <= N_RUNS; irun++) {
        pi<<<GRID_SIZE, BLOCK_SIZE>>>(d_pis, states, xs, radii, N_TRIALS);
    }
    gpuErrchk(cudaMemcpy(h_pis, d_pis, N_KERNELS*sizeof(double), cudaMemcpyDeviceToHost));
    double avg = 0.0;
    for(int i=0; i < N_KERNELS; i++) {
        avg += h_pis[i] / N_RUNS;
    }
    avg /= N_KERNELS;

    time_t end = clock();

    int64_t iters = int64_t(N_KERNELS)*int64_t(N_TRIALS)*int64_t(N_RUNS);
    int elapsed = 1000*(end-start)/CLOCKS_PER_SEC;

    cout.precision(10);
    cout << "pi = " << 2*avg << "\n";
    cout << elapsed << " ms" << "\n";
    cout << double(iters)/double(elapsed) << " iters/ms\n";

    free(h_pis);
    cudaFree(d_pis);
    cudaFree(states);
    cudaFree(xs);
    return 0;
}
