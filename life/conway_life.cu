#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "conway_life.hpp"

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

void freeGameBoard(GameBoard * g) {
    free(g->board);
    free(g->work);
    free(g);
}


char cellCharacter(uint8_t status) {
    return (status == ALIVE)? '#' : '.';
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
    int i = g->R/2;
    int j = g->C/2;
    g->board[g->C*(i+0)+(j+1)] = ALIVE;
    g->board[g->C*(i+0)+(j+2)] = ALIVE;
    g->board[g->C*(i+1)+(j+0)] = ALIVE;
    g->board[g->C*(i+1)+(j+1)] = ALIVE;
    g->board[g->C*(i+2)+(j+1)] = ALIVE;

    time_t start = clock();
    updateBoard(g, t);
    time_t end = clock();
    int elapsed = 1000*(end-start)/CLOCKS_PER_SEC;

    printBoard(g);
    printf("Took %d ms\n", elapsed);

    freeGameBoard(g);
    return 0;
}
