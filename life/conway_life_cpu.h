#ifndef CONWAY_LIFE_CPU_H
#define CONWAY_LIFE_CPU_H

#include <stdint.h>

uint8_t const DEAD  = 0;
uint8_t const ALIVE = 1;

typedef struct {
    uint8_t * board;
    uint8_t * work;
    size_t R, C;
} GameBoard;

GameBoard * newGameBoard(size_t R, size_t C);
void printBoard(GameBoard * g);
uint8_t newStatus(uint8_t status, int count);
char cellCharacter(uint8_t status);
void updateCell(uint8_t * next, uint8_t * cur, int i, int j, size_t R, size_t C);
void updateBoard(GameBoard * g);
void displayBoard(GameBoard * g);
void printBoard(GameBoard * g);
void runGame(GameBoard * g);

#endif
