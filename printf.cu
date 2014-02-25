// printf.cu
/*
 * Simple script to show how to print on a GPU.
 * NOTE: I have no idea why, but this simply does not work unless
 * you initialize an array on the GPU, then free it afterwards.
 * Hence the cudaMalloc() and cudaFree() calls that do nothing.
 *
 * My assumption is that it has something to do with initializing
 * an "active" state on the device. cudaDeviceSynchronize(),
 * cudaDeviceReset(), and cudaSetDevice() do not suffice, however.
 */

#include <stdio.h>

__global__ void hello() {
    printf("Hello from Block %d, Thread %d\n", blockIdx.x, threadIdx.x);
}

int main() {
    float *d_arr;
    cudaMalloc(&d_arr, 25*sizeof(float));

    hello<<<5,5>>>();

    cudaFree(d_arr);
    return 0;
}
