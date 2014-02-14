// reduce.cu
/*
 * An exploration of various reduce algorithms on a CUDA GPU
 */

#include <stdio.h>
#include <assert.h>
#include <cuda_runtime.h>

// Straightforward parallel reduce
__global__ void reduce1(float* d_out, float* d_in, int n) {
    // sizeof(d_out) == gridDim.x ((one spot for the sum of each block))
    int t_idx = threadIdx.x;
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    extern __shared__ float sdata[]; // dynamically allocated shared memory

    sdata[t_idx] = d_in[idx];
    __syncthreads();

    // Start by having every even# element add its neighbor
    // Then renumber everything and do it again (the s *= 2 is "renumbering")
    for(int s=1; s < blockDim.x; s *= 2) {
        if( t_idx % (2*s) == 0 ) {
            sdata[t_idx] += sdata[t_idx + s];
        }
        __syncthreads();
    }

    // At the end, the first thread in every block will have the sum of that block
    if( t_idx == 0 ) {
        d_out[blockIdx.x] = sdata[0];
    }
}

// Same as reduce1(), but with a different summation scheme
__global__ void reduce2(float* d_out, float* d_in, int n) {
    // sizeof(d_out) == gridDim.x ((one spot for the sum of each block))
    int t_idx = threadIdx.x;
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    extern __shared__ float sdata[]; // dynamically allocated shared memory

    sdata[t_idx] = d_in[idx];
    __syncthreads();

    // Start by folding the second half of the block into the first half
    // Then fold the second quarter into the first quarter
    // Then the second eighth into the first, etc
    assert( (blockDim.x & (blockDim.x-1)) == 0 ); // for now, make sure it's a power of 2
    for(int s=blockDim.x/2; s > 0; s /= 2) {
        if( t_idx < s ) { // if t_idx is in the first {half/quarter/etc}
            sdata[t_idx] += sdata[t_idx + s];
        }
        __syncthreads();
    }
    if( t_idx == 0 ) {
        d_out[blockIdx.x] = sdata[0];
    }
}

// Same as reduce2(), but each block takes care of chunks of the input array
// evenly spaced out at grimDim.x*blockDim.x spaces
__global__ void reduce3(float* d_out, float* d_in, int n) {
    // sizeof(d_out) == gridDim.x ((one spot for the sum of each block))
    int t_idx = threadIdx.x;
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    extern __shared__ float sdata[]; // dynamically allocated shared memory

    sdata[t_idx] = 0;
    for(int i=idx; i < n; i += gridDim.x*blockDim.x) {
        sdata[t_idx] += d_in[i];
    }
    __syncthreads();

    // Start by folding the second half of the block into the first half
    // Then fold the second quarter into the first quarter
    // Then the second eighth into the first, etc
    assert( (blockDim.x & (blockDim.x-1)) == 0 ); // for now, make sure it's a power of 2
    for(int s=blockDim.x/2; s > 0; s /= 2) {
        if( t_idx < s ) { // if t_idx is in the first {half/quarter/etc}
            sdata[t_idx] += sdata[t_idx + s];
        }
        __syncthreads();
    }
    if( t_idx == 0 ) {
        d_out[blockIdx.x] = sdata[0];
    }
}


double testReduce(int n) {
    const int BLOCK_SIZE = 512;
    const int NUM_BLOCKS = n/BLOCK_SIZE;
    assert(n%BLOCK_SIZE == 0);

    float *h_in;
    float *h_out;
    h_in  = (float*)malloc(n*sizeof(float));
    h_out = (float*)malloc(NUM_BLOCKS*sizeof(float));

    double ans = 0.0;
    for(int i=0; i < n; i++) {
        h_in[i] = float((i%3) - 1);
    }

    float *d_in;
    float *d_out;
    cudaMalloc(&d_in, n*sizeof(float));
    cudaMalloc(&d_out, NUM_BLOCKS*sizeof(float));

    cudaMemcpy(d_in, h_in, n*sizeof(float), cudaMemcpyHostToDevice);

    clock_t start = clock();
    reduce1<<<NUM_BLOCKS,BLOCK_SIZE,BLOCK_SIZE*sizeof(float)>>>(d_out, d_in, n);
    clock_t end = clock();

    cudaMemcpy(h_out, d_out, NUM_BLOCKS*sizeof(float), cudaMemcpyDeviceToHost);

    double sum = 0.0;
    for(int i=0; i < NUM_BLOCKS; i++) {
        sum += h_out[i];
    }

    cudaFree(d_in);
    cudaFree(d_out);
    free(h_in);
    free(h_out);

    return (end-start)/double(CLOCKS_PER_SEC);
}

int main() {
    double tot = 0.0;
    for(int i=0; i < 500; i++) {
        tot += testReduce(1<<18);
    }
    printf("tot = %g\n", tot);
    return 0;
}
