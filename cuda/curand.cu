#include <stdio.h>
#include <curand.h>

int main() {
    int n = 20;
    float* h_xs;
    float* d_xs;

    h_xs = (float*)malloc(n*sizeof(float));
    cudaMalloc(&d_xs, n*sizeof(float));


    curandGenerator_t prng;
    curandCreateGenerator(&prng, CURAND_RNG_PSEUDO_MTGP32); // single-precision
    curandSetPseudoRandomGeneratorSeed(prng, 42);
    curandGenerateUniform(prng, d_xs, n);
    curandDestroyGenerator(prng);

    cudaMemcpy(h_xs, d_xs, n*sizeof(float), cudaMemcpyDeviceToHost);
    for(int i=0; i < n; i++) {
        printf("%f\n", h_xs[i]);
    }

    cudaFree(d_xs);
    free(h_xs);

    return 0;
}
