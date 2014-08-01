#include <stdio.h>

int main() {
    int count;
    cudaGetDeviceCount(&count);

    printf("%d\n", count);

    return 0;
}
