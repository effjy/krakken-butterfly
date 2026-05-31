#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "krakken.h"

// Define krakken_stream_t and the parallel function prototype
typedef struct {
  uint8_t *out;
  size_t outlen;
  const uint8_t *in;
  size_t inlen;
} krakken_stream_t;

int krakken_hash_avx2_parallel(const krakken_stream_t *streams, int n, int n_threads);

// Pseudorandom data generator
static void fill_random(uint8_t *buf, size_t len) {
    uint32_t state = 0x12345678;
    for (size_t i = 0; i < len; i++) {
        state ^= state << 13;
        state ^= state >> 17;
        state ^= state << 5;
        buf[i] = (uint8_t)(state & 0xFF);
    }
}

static void print_hash(const uint8_t *hash, size_t len) {
    for (size_t i = 0; i < len; i++) {
        printf("%02x", hash[i]);
    }
    printf("\n");
}

int main(void) {
    printf("=== Starting Krakken Verification Tests ===\n");

    size_t lengths[] = {0, 1, 5, 10, 158, 159, 160, 161, 319, 320, 321, 1000, 4096};
    int num_lengths = sizeof(lengths) / sizeof(lengths[0]);
    size_t out_lengths[] = {16, 32, 64};
    int num_out_lengths = sizeof(out_lengths) / sizeof(out_lengths[0]);

    // Test single-stream scalar vs AVX2
    for (int i = 0; i < num_lengths; i++) {
        size_t inlen = lengths[i];
        uint8_t *in = malloc(inlen > 0 ? inlen : 1);
        if (inlen > 0) {
            fill_random(in, inlen);
        }

        for (int j = 0; j < num_out_lengths; j++) {
            size_t outlen = out_lengths[j];
            uint8_t *out_scalar = malloc(outlen);
            uint8_t *out_avx2 = malloc(outlen);

            memset(out_scalar, 0, outlen);
            memset(out_avx2, 0, outlen);

            krakken_hash_scalar(out_scalar, outlen, in, inlen);
            krakken_hash_avx2(out_avx2, outlen, in, inlen);

            if (memcmp(out_scalar, out_avx2, outlen) != 0) {
                printf("FAIL: Mismatch at input length %zu, output length %zu\n", inlen, outlen);
                printf("Scalar: ");
                print_hash(out_scalar, outlen);
                printf("AVX2  : ");
                print_hash(out_avx2, outlen);
                free(in);
                free(out_scalar);
                free(out_avx2);
                exit(1);
            }

            free(out_scalar);
            free(out_avx2);
        }
        free(in);
    }
    printf("SUCCESS: All single-stream tests matched perfectly!\n");

    // Test parallel implementation
    printf("Testing parallel AVX2 hashing...\n");
    int num_streams = 32;
    krakken_stream_t *streams = malloc(num_streams * sizeof(krakken_stream_t));
    uint8_t **inputs = malloc(num_streams * sizeof(uint8_t*));
    uint8_t **outputs_parallel = malloc(num_streams * sizeof(uint8_t*));
    uint8_t **outputs_scalar = malloc(num_streams * sizeof(uint8_t*));

    for (int i = 0; i < num_streams; i++) {
        size_t inlen = (i * 97) % 2000; // varying input sizes
        size_t outlen = 32;
        inputs[i] = malloc(inlen > 0 ? inlen : 1);
        if (inlen > 0) {
            fill_random(inputs[i], inlen);
        }
        outputs_parallel[i] = malloc(outlen);
        outputs_scalar[i] = malloc(outlen);

        streams[i].in = inputs[i];
        streams[i].inlen = inlen;
        streams[i].out = outputs_parallel[i];
        streams[i].outlen = outlen;

        krakken_hash_scalar(outputs_scalar[i], outlen, inputs[i], inlen);
    }

    int ret = krakken_hash_avx2_parallel(streams, num_streams, 4);
    assert(ret == 0);

    for (int i = 0; i < num_streams; i++) {
        if (memcmp(outputs_scalar[i], outputs_parallel[i], 32) != 0) {
            printf("FAIL: Parallel mismatch at stream index %d\n", i);
            exit(1);
        }
    }
    printf("SUCCESS: Parallel AVX2 tests matched scalar results perfectly!\n");

    // Cleanup parallel test memory
    for (int i = 0; i < num_streams; i++) {
        free(inputs[i]);
        free(outputs_parallel[i]);
        free(outputs_scalar[i]);
    }
    free(streams);
    free(inputs);
    free(outputs_parallel);
    free(outputs_scalar);

    printf("=== Krakken Verification Completed Successfully ===\n");
    return 0;
}
