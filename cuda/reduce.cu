// reduce.cu
/*
 * An exploration of various reduce algorithms on a CUDA GPU
 */

#include <stdio.h>
#include <assert.h>
#include <cuda_runtime.h>

__global__ void reduce_atomic(float* d_arr, int len) {
    for(int idx=blockIdx.x*blockDim.x + threadIdx.x; idx < len; idx += gridDim.x*blockDim.x) {
        if(idx == 0) { continue; }
        atomicAdd(&d_arr[0], d_arr[idx]);
    }
}



int main() {
    int n = 1024;
    int bytes = n * sizeof(float);
    float ans = n*(n+1)/2;

    float *h_arr;
    h_arr = (float*)malloc(bytes);
    for(int i=0; i < n; i++) {
        h_arr[i] = float(i+1);
    }

    float *d_arr;
    cudaMalloc(&d_arr, bytes);
    cudaMemcpy(d_arr, h_arr, bytes, cudaMemcpyHostToDevice);

    reduce_atomic<<<5,5>>>(d_arr, n);

    cudaMemcpy(h_arr, d_arr, sizeof(float), cudaMemcpyDeviceToHost);
    if( h_arr[0] == ans ) {
        printf("[PASS]");
    } else {
        printf("[!!!FAIL!!!]");
    }
    printf(" sum = %f\n", h_arr[0]);

    cudaFree(d_arr);
    free(h_arr);
    return 0;
}
