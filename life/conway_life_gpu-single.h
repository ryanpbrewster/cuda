#ifndef CONWAY_LIFE_GPU_SINGLE_H
#define CONWAY_LIFE_GPU_SINGLE_H

#include <stdint.h>

uint8_t const DEAD  = 0;
uint8_t const ALIVE = 1;

typedef struct {
    uint8_t * board;
    uint8_t * work;
    size_t R, C;
} GameBoard;

__device__ uint8_t newStatus(uint8_t status, int count);
__global__ void updateCell(uint8_t * next, uint8_t * cur, size_t R, size_t C);

__host__ GameBoard * newGameBoard(size_t R, size_t C);
__host__ void printBoard(GameBoard * g);
__host__ char cellCharacter(uint8_t status);
__host__ void updateBoard(GameBoard * g);
__host__ void displayBoard(GameBoard * g);
__host__ void printBoard(GameBoard * g);
__host__ void runGame(GameBoard * g);

#endif
