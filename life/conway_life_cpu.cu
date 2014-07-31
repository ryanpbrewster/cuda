#include "conway_life.hpp"

void updateCell(uint8_t * next, uint8_t * cur, int i, int j, size_t R, size_t C) {
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

void updateBoard(GameBoard * g, int t) {
    for(int gen=1; gen <= t; gen++) {
        for(int i=0; i < g->R; i++) {
            for(int j=0; j < g->C; j++) {
                updateCell(g->work, g->board, i, j, g->R, g->C);
            }
        }

        uint8_t * tmp = g->board;
        g->board = g->work;
        g->work = tmp;
    }
}
