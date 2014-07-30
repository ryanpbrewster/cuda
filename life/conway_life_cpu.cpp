#include <iostream>
#include <stdlib.h>
#include <string.h>

using namespace std;

const uint8_t DEAD  = 0;
const uint8_t ALIVE = 1;

uint8_t updateCell(uint8_t cur, int count) {
    if( cur == ALIVE ) {
        return (count < 2 || count > 3)? DEAD : ALIVE;
    } else {
        return (count == 3)? ALIVE : DEAD;
    }
}

char cellCharacter(uint8_t val) {
    return (val == ALIVE)? '#' : '.';
}

void updateBoard(uint8_t ** boardptr, uint8_t ** workptr, size_t const R, size_t const C) {
    uint8_t * board = *boardptr;
    uint8_t * work = *workptr;

    for(int i=0; i < R; i++) {
        for(int j=0; j < R; j++) {
            // neighbors of (i,j)
            int ns[] = { C*((i-1)%R)  + ((j-1)%C)
                       , C*((i-1)%R)  + ((j  )%C)
                       , C*((i-1)%R)  + ((j+1)%C)
                       , C*((i  )%R)  + ((j-1)%C)
                       , C*((i  )%R)  + ((j+1)%C)
                       , C*((i+1)%R)  + ((j-1)%C)
                       , C*((i+1)%R)  + ((j  )%C)
                       , C*((i+1)%R)  + ((j+1)%C)
                       };
            int idx = C*i + j;
            int count = board[ns[0]] + board[ns[1]] + board[ns[2]] + board[ns[3]] + board[ns[4]] + board[ns[5]] + board[ns[6]] + board[ns[7]];
            work[idx] = updateCell(board[idx], count);
        }
    }

    uint8_t * tmp = *boardptr;
    *boardptr = *workptr;
    *workptr = tmp;
}

void printBoard(uint8_t const * board, size_t const R, size_t const C) {
    for(int i=0; i < R; i++) {
        for(int j=0; j < C; j++) {
            cout << cellCharacter(board[C*i+j]);
        }
        cout << "\n";
    }
}

int main() {
    size_t const R = 10;
    size_t const C = 20;
    size_t const BOARD_BYTES = R*C*sizeof(uint8_t);
    uint8_t* board = (uint8_t*) malloc(BOARD_BYTES);
    uint8_t* work  = (uint8_t*) malloc(BOARD_BYTES);
    memset(board, DEAD, BOARD_BYTES);

    board[C*0+1] = ALIVE;
    board[C*1+1] = ALIVE;
    board[C*2+1] = ALIVE;

    printBoard(board, R, C);
    cout << "\n\n";
    updateBoard(&board, &work, R, C);
    printBoard(board, R, C);
    cout << "\n\n";
    updateBoard(&board, &work, R, C);
    printBoard(board, R, C);

    return 0;
}
