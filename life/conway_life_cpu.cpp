#include <ncurses.h>
#include <iostream>
#include <stdlib.h>
#include <string.h>

using namespace std;

const uint8_t DEAD  = 0;
const uint8_t ALIVE = 1;

typedef struct {
    uint8_t * board;
    uint8_t * work;
    size_t R, C;
} GameBoard;

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

void printBoard(GameBoard * g);

uint8_t newStatus(uint8_t status, int count) {
    if( status == ALIVE ) {
        return (count < 2 || count > 3)? DEAD : ALIVE;
    } else {
        return (count == 3)? ALIVE : DEAD;
    }
}

char cellCharacter(uint8_t status) {
    return (status == ALIVE)? '#' : '.';
}

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

void updateBoard(GameBoard * g) {
    for(int i=0; i < g->R; i++) {
        for(int j=0; j < g->C; j++) {
            updateCell(g->work, g->board, i, j, g->R, g->C);
        }
    }

    uint8_t * tmp = g->board;
    g->board = g->work;
    g->work = tmp;
}

void displayBoard(GameBoard * g) {
    for(int i=0; i < g->R; i++) {
        for(int j=0; j < g->C; j++) {
            addch(cellCharacter(g->board[g->C*i+j]));
        }
        addch('\n');
    }
}
void printBoard(GameBoard * g) {
    for(int i=0; i < g->R; i++) {
        for(int j=0; j < g->C; j++) {
            cout << cellCharacter(g->board[g->C*i+j]);
        }
        cout << "\n";
    }
}


void runGame(GameBoard * g) {
    displayBoard(g);
    refresh();
    while(true) {
        char ch = getch();
        switch(ch) {
            case 'q': return;
            default:  updateBoard(g); break;
        }
        clear();
        displayBoard(g);
        refresh();
    }
}

int main() {
    initscr();
    cbreak();
    noecho();

    GameBoard * g = newGameBoard(15, 60);

    // Create the R-pentomino somewhere near the center of the board
    int i = g->R/2;
    int j = g->C/2;
    g->board[g->C*(i+0)+(j+1)] = ALIVE;
    g->board[g->C*(i+0)+(j+2)] = ALIVE;
    g->board[g->C*(i+1)+(j+0)] = ALIVE;
    g->board[g->C*(i+1)+(j+1)] = ALIVE;
    g->board[g->C*(i+2)+(j+1)] = ALIVE;

    runGame(g);

    endwin();

    free(g->board);
    free(g->work);
    free(g);
    return 0;
}
