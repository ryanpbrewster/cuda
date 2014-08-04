// metropolis1.cu
/*
 * A simple CUDA-enabled program that approximates Pi by evaluating
 *     Integrate[ Sqrt[1-x^2], {x,-1,1} ]
 * using Metropolis Monte Carlo, with the weight function A*(1-x^2)
 * where A is a normalization factor
 */

#include <iostream>
#include <curand.h>
#include <curand_kernel.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdbool.h>
#include "metropolis.hpp"
using namespace std;

int main(int argc, char** argv) {
    int GRID_SIZE  = 256;
    int BLOCK_SIZE = 256;
    int N_TRIALS   = 1000;
    int N_RUNS     = 10;

    char x;
    opterr = 0;
    while((x = getopt(argc, argv, "g:b:t:r:")) != -1) {
        switch(x) {
            case 'g': GRID_SIZE  = atoi(optarg); break;
            case 'b': BLOCK_SIZE = atoi(optarg); break;
            case 't': N_TRIALS   = atoi(optarg); break;
            case 'r': N_RUNS     = atoi(optarg); break;
            case '?':
                if (optopt == 'g' || optopt == 'b' || optopt == 't' || optopt == 'r') {
                    fprintf (stderr, "Option -%c requires an argument.\n", optopt);
                } else {
                    fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                }
                abort();
            default: abort();
        }
    }

    printf("G=% 4d, B=% 4d, T=% 8d, R=% 5d\n", GRID_SIZE, BLOCK_SIZE, N_TRIALS, N_RUNS);

    MPI_Init(NULL, NULL);
    int NPROC;
    MPI_Comm_size(MPI_COMM_WORLD, &NPROC);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double local_pi = mcpi(GRID_SIZE, BLOCK_SIZE, N_TRIALS, N_RUNS, GRID_SIZE*BLOCK_SIZE*rank);
    double global_pi;
    MPI_Reduce(&local_pi, &global_pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Finalize();

    printf("%f\n", global_pi);
}
