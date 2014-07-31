#include "conway_life.hpp"

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, char *file, int line, bool abort=true) {
    if (code != cudaSuccess) {
        fprintf(stderr,"GPUassert: \"%s\" in %s at %d\n", cudaGetErrorString(code), file, line);
        if (abort) {
            exit(code);
        }
    }
}

__global__ void updateCell(uint8_t * next, uint8_t * cur, size_t R, size_t C) {
    for(int i=blockDim.y*blockIdx.y + threadIdx.y; i < R; i += blockDim.y*gridDim.y) {
        for(int j=blockDim.x*blockIdx.x + threadIdx.x; j < C; j += blockDim.x*gridDim.x) {
            // neighbors of (i,j)
            int ns[] = { C*((i-1+R)%R)  + ((j-1+C)%C)
                       , C*((i-1+R)%R)  + ((j  +C)%C)
                       , C*((i-1+R)%R)  + ((j+1+C)%C)
                       , C*((i  +R)%R)  + ((j-1+C)%C)
                       , C*((i  +R)%R)  + ((j+1+C)%C)
                       , C*((i+1+R)%R)  + ((j-1+C)%C)
                       , C*((i+1+R)%R)  + ((j  +C)%C)
                       , C*((i+1+R)%R)  + ((j+1+C)%C)
                       };
            int idx = C*i + j;
            int count = cur[ns[0]] + cur[ns[1]] + cur[ns[2]] + cur[ns[3]] + cur[ns[4]] + cur[ns[5]] + cur[ns[6]] + cur[ns[7]];
            next[idx] = newStatus(cur[idx], count);
        }
    }
}

void updateBoard(GameBoard * g, int t) {
    uint8_t * d_board;
    uint8_t * d_work;
    size_t const bytes = g->R * g->C * sizeof(uint8_t);
    gpuErrchk( cudaMalloc(&d_board, bytes) );
    gpuErrchk( cudaMalloc(&d_work,  bytes) );

    gpuErrchk( cudaMemcpy(d_board, g->board, bytes, cudaMemcpyHostToDevice) );

    for(int gen=1; gen <= t; gen++) {
        updateCell<<<dim3(16,16,1), dim3(16,16,1)>>>(d_work, d_board, g->R, g->C);
        uint8_t * tmp = d_board;
        d_board = d_work;
        d_work = tmp;
    }
    gpuErrchk( cudaMemcpy(g->board, d_board, bytes, cudaMemcpyDeviceToHost) );
    gpuErrchk( cudaFree(d_board) );
    gpuErrchk( cudaFree(d_work) );
}
