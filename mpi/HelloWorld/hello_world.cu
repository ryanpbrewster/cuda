#include "hello_world.hpp"

int countGPUs() {
    int gpu_count;
    cudaGetDeviceCount(&gpu_count);
    return gpu_count;
}
