#ifndef CONWAY_LIFE_HPP
#define CONWAY_LIFE_HPP

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

static uint8_t const DEAD  = 0;
static uint8_t const ALIVE = 1;

typedef struct {
    uint8_t * board;
    uint8_t * work;
    size_t R, C;
} GameBoard;

__host__ __device__
inline uint8_t newStatus(uint8_t status, int count) {
    if( status == ALIVE ) {
        return (count < 2 || count > 3)? DEAD : ALIVE;
    } else {
        return (count == 3)? ALIVE : DEAD;
    }
}

void updateBoard(GameBoard * g, int t);

#endif
