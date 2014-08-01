#include <mpi.h>
#include <stdio.h>
#include "hello_world.hpp"

int main(int argc, char** argv) {
  // Initialize the MPI environment
  MPI_Init(NULL, NULL);

  // Get the number of processes
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // Get the rank of the process
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  int gpu_count = countGPUs();

  // Print off a hello world message
  printf("Processor %d of %d checking in with %d gpus\n", world_rank, world_size, gpu_count);

  // Finalize the MPI environment.
  MPI_Finalize();
}
