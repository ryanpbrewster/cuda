#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "conway_life_gpu-single.h"

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, char *file, int line, bool abort=true) {
    if (code != cudaSuccess) {
        fprintf(stderr,"GPUassert: \"%s\" in %s at %d\n", cudaGetErrorString(code), file, line);
        if (abort) {
            exit(code);
        }
    }
}

void freeGameBoard(GameBoard * g) {
    free(g->board);
    free(g->work);
    free(g);
}

GameBoard * newGameBoard(size_t R, size_t C) {
    GameBoard * g = (GameBoard *) malloc(sizeof(GameBoard));
    g->R = R;
    g->C = C;
    size_t const BOARD_BYTES = R*C*sizeof(uint8_t);
    g->board = (uint8_t*) malloc(BOARD_BYTES);
    g->work  = (uint8_t*) malloc(BOARD_BYTES);
    memset(g->board, DEAD, BOARD_BYTES);
    return g;
}

__device__ uint8_t newStatus(uint8_t status, int count) {
    if( status == ALIVE ) {
        return (count < 2 || count > 3)? DEAD : ALIVE;
    } else {
        return (count == 3)? ALIVE : DEAD;
    }
}

char cellCharacter(uint8_t status) {
    return (status == ALIVE)? '#' : '.';
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

void updateBoard(GameBoard * g, uint8_t * d_board, uint8_t * d_work) {
    size_t const bytes = g->R * g->C * sizeof(uint8_t);
    gpuErrchk( cudaMemcpy(d_board, g->board, bytes, cudaMemcpyHostToDevice) );
    updateCell<<<dim3(16,16,1), dim3(16,16,1)>>>(d_work, d_board, g->R, g->C);
    gpuErrchk( cudaMemcpy(g->board, d_work, bytes, cudaMemcpyDeviceToHost) );
}

void printBoard(GameBoard * g) {
    for(int i=0; i < g->R; i++) {
        for(int j=0; j < g->C; j++) {
            printf("%c", cellCharacter(g->board[g->C*i+j]));
        }
        printf("\n");
    }
}

int main(int argc, char** argv) {
    if( argc != 2 ) {
        fprintf(stderr, "Usage: %s [# of generations to simulate]\n", argv[0]);
        return 1;
    }
    int t = atoi(argv[1]);

    GameBoard * g = newGameBoard(15, 60);

    // Create the R-pentomino somewhere near the center of the board
    int mi = g->R/2;
    int mj = g->C/2;
    g->board[g->C*(mi+0)+(mj+1)] = ALIVE;
    g->board[g->C*(mi+0)+(mj+2)] = ALIVE;
    g->board[g->C*(mi+1)+(mj+0)] = ALIVE;
    g->board[g->C*(mi+1)+(mj+1)] = ALIVE;
    g->board[g->C*(mi+2)+(mj+1)] = ALIVE;

    uint8_t * d_board;
    uint8_t * d_work;
    size_t const bytes = g->R * g->C * sizeof(uint8_t);
    gpuErrchk( cudaMalloc(&d_board, bytes) );
    gpuErrchk( cudaMalloc(&d_work,  bytes) );

    for(int i=1; i <= t; i++) {
        updateBoard(g, d_board, d_work);
    }
    printBoard(g);

    gpuErrchk( cudaFree(d_board) );
    gpuErrchk( cudaFree(d_work) );
    freeGameBoard(g);
    return 0;
}
