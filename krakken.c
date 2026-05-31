#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <pthread.h>
#include "krakken.h"

#define KRAKKEN_ROUNDS 8

static void krakken_memset_wrap(void *p, int c, size_t n) { memset(p, c, n); }
static void (* volatile krakken_memset_volatile)(void *, int, size_t) = krakken_memset_wrap;

static void krakken_secure_zero(void *ptr, size_t n) {
    krakken_memset_volatile(ptr, 0, n);
}

static inline uint64_t rotl64(uint64_t x, int n) {
    n &= 63;
    if (n == 0) return x;
    return (x << n) | (x >> (64 - n));
}

static inline uint64_t rotr64(uint64_t x, int n) {
    n &= 63;
    if (n == 0) return x;
    return (x >> n) | (x << (64 - n));
}

static inline uint64_t custom_sbox8_64(uint64_t w) {
    uint64_t res = 0;
    for (int i = 0; i < 256; i++) {
        uint64_t entry = ABYSSAL_SBOX[i];
        uint64_t val_i = i * 0x0101010101010101ULL;
        uint64_t diff = w ^ val_i;
        uint64_t temp = (diff & 0x7F7F7F7F7F7F7F7FULL) + 0x7F7F7F7F7F7F7F7FULL;
        uint64_t zero_msb = ~(temp | diff) & 0x8080808080808080ULL;
        uint64_t mask = zero_msb >> 7;
        res |= mask * entry;
    }
    return res;
}

static const int rho[32] = {
    32,  1, 62, 28, 36, 44, 15, 61,
     6, 19, 24, 55,  3, 10, 43, 17,
    25, 39, 41, 59, 47,  8, 56, 14,
    18, 35, 21, 33,  2, 49, 22, 51
};

#define KECCAK_ROUNDS 24
static const uint64_t keccak_rc[24] = {
    0x0000000000000001ULL, 0x0000000000008082ULL,
    0x800000000000808AULL, 0x8000000080008000ULL,
    0x000000000000808BULL, 0x0000000080000001ULL,
    0x8000000080008081ULL, 0x8000000000008009ULL,
    0x000000000000008AULL, 0x0000000000000088ULL,
    0x0000000080008009ULL, 0x000000008000000AULL,
    0x000000008000808BULL, 0x800000000000008BULL,
    0x8000000000008089ULL, 0x8000000000008003ULL,
    0x8000000000008002ULL, 0x8000000000000080ULL,
    0x000000000000800AULL, 0x800000008000000AULL,
    0x8000000080008081ULL, 0x8000000000008080ULL,
    0x0000000080000001ULL, 0x8000000080008008ULL
};

static const int keccak_rho[24] = {
     1,  3,  6, 10, 15, 21, 28, 36, 45, 55,  2, 14,
    27, 41, 56,  8, 25, 43, 62, 18, 39, 61, 20, 44
};

static const int keccak_pi[24] = {
    10,  7, 11, 17, 18,  3,  5, 16,  8, 21, 24,  4,
    15, 23, 19, 13, 12,  2, 20, 14, 22,  9,  6,  1
};

static void keccakf1600(uint64_t st[25]) {
    for (int r = 0; r < KECCAK_ROUNDS; r++) {
        uint64_t bc[5], t;
        for (int i = 0; i < 5; i++)
            bc[i] = st[i] ^ st[i+5] ^ st[i+10] ^ st[i+15] ^ st[i+20];
        for (int i = 0; i < 5; i++) {
            t = bc[(i+4)%5] ^ rotl64(bc[(i+1)%5], 1);
            for (int j = 0; j < 25; j += 5) st[j+i] ^= t;
        }
        t = st[1];
        for (int i = 0; i < 24; i++) {
            int j = keccak_pi[i];
            bc[0] = st[j];
            st[j] = rotl64(t, keccak_rho[i]);
            t = bc[0];
        }
        for (int j = 0; j < 25; j += 5) {
            for (int i = 0; i < 5; i++) bc[i] = st[j+i];
            for (int i = 0; i < 5; i++)
                st[j+i] ^= (~bc[(i+1)%5]) & bc[(i+2)%5];
        }
        st[0] ^= keccak_rc[r];
    }
}

static void shake128_squeeze(const char *domain, uint8_t *out, size_t outlen) {
    uint64_t st[25] = {0};
    size_t rate = 168;
    uint8_t *st8 = (uint8_t *)st;

    size_t dlen = strlen(domain);
    if (dlen >= rate - 2) {
        fprintf(stderr, "FATAL: Domain string too long in shake128_squeeze\n");
        abort();
    }
    for (size_t i = 0; i < dlen; i++) st8[i] ^= (uint8_t)domain[i];
    st8[dlen]     ^= 0x1F;
    st8[rate - 1] ^= 0x80;
    keccakf1600(st);

    size_t done = 0;
    while (done < outlen) {
        size_t take = outlen - done < rate ? outlen - done : rate;
        memcpy(out + done, st8, take);
        done += take;
        if (done < outlen) keccakf1600(st);
    }
}

static uint64_t rc[KRAKKEN_ROUNDS][32];
static pthread_once_t rc_once = PTHREAD_ONCE_INIT;

static void _rc_init_impl(void) {
    uint8_t buf[KRAKKEN_ROUNDS * 32 * 8];
    shake128_squeeze("Krakken-2048 Abyssal v1 - Primary ", buf, sizeof(buf));
    for (int ir = 0; ir < KRAKKEN_ROUNDS; ir++) {
        for (int i = 0; i < 32; i++) {
            const uint8_t *p = buf + (ir * 32 + i) * 8;
            uint64_t v = (uint64_t)p[0]
                       | ((uint64_t)p[1] <<  8)
                       | ((uint64_t)p[2] << 16)
                       | ((uint64_t)p[3] << 24)
                       | ((uint64_t)p[4] << 32)
                       | ((uint64_t)p[5] << 40)
                       | ((uint64_t)p[6] << 48)
                       | ((uint64_t)p[7] << 56);
            rc[ir][i] = v ? v : 0xDEADBEEFCAFEBABEULL;
        }
    }
}

void init_rc_vectors(void) {
    pthread_once(&rc_once, _rc_init_impl);
}

const uint64_t *rc_get(int round) { return rc[round]; }

void theta_scalar(uint64_t state[32]) {
    uint64_t parity[8];
    for (int c = 0; c < 8; c++)
        parity[c] = state[4*c] ^ state[4*c+1] ^ state[4*c+2] ^ state[4*c+3];
    for (int c = 0; c < 8; c++) {
        uint64_t d = rotr64(parity[(c + 7) & 7], 1) ^ parity[(c + 1) & 7];
        for (int y = 0; y < 4; y++)
            state[4*c + y] ^= d;
    }
}

static inline uint64_t gf28_double_word(uint64_t w) {
    uint64_t shifted = (w << 1) & 0xFEFEFEFEFEFEFEFEULL;
    uint64_t msb_mask = (w & 0x8080808080808080ULL) >> 7;
    uint64_t reduction = msb_mask * 0x1DULL;
    return shifted ^ reduction;
}

static inline uint64_t gf28_mul_word_fast(uint64_t w, uint8_t k) {
    switch (k) {
        case 0x01: return w;
        case 0x02: return gf28_double_word(w);
        case 0x04: return gf28_double_word(gf28_double_word(w));
        case 0x05: {
            uint64_t w2 = gf28_double_word(w);
            return gf28_double_word(w2) ^ w;
        }
        case 0x08: return gf28_double_word(gf28_double_word(gf28_double_word(w)));
        case 0x09: {
            uint64_t w2 = gf28_double_word(w);
            uint64_t w4 = gf28_double_word(w2);
            return gf28_double_word(w4) ^ w;
        }
        default:
            return w;
    }
}

static const uint8_t mds_coeffs[8] = { 0x01, 0x01, 0x04, 0x01, 0x08, 0x05, 0x02, 0x09 };

void tentacle_mds_scalar(uint64_t state[32]) {
    for (int y = 0; y < 4; y++) {
        uint64_t row[8];
        for (int c = 0; c < 8; c++) row[c] = state[c*4 + y];
        for (int c = 0; c < 8; c++) {
            uint64_t sum = 0;
            for (int i = 0; i < 8; i++) {
                sum ^= gf28_mul_word_fast(row[(c + i) & 7], mds_coeffs[i]);
            }
            state[c*4 + y] = sum;
        }
    }
}

void rho_scalar(uint64_t state[32]) {
    for (int i = 0; i < 32; i++)
        state[i] = rotl64(state[i], rho[i]);
}

void pi_scalar(uint64_t state[32]) {
    uint64_t temp[32];
    for (int i = 0; i < 32; i++) {
        int x = i / 4, y = i % 4;
        int new_x = (x + 3 * y) & 7;
        temp[new_x * 4 + y] = state[i];
    }
    memcpy(state, temp, 32 * sizeof(uint64_t));
}

void chi_scalar(uint64_t state[32]) {
    for (int y = 0; y < 4; y++) {
        for (int p = 0; p < 4; p++) {
            int ca = (p * 2)     * 4 + y;
            int cb = (p * 2 + 1) * 4 + y;
            uint64_t a = state[ca], b = state[cb];
            uint64_t ap = custom_sbox8_64(a ^ rotl64(b,  32));
            uint64_t bp = custom_sbox8_64(b ^ rotl64(ap, 32));
            state[ca] = ap;
            state[cb] = bp;
        }
    }
}

void butterfly_diffusion_scalar(uint64_t state[32]) {
    static const int rotations[5] = { 13, 23, 37, 41, 53 };
    for (int stage = 0; stage < 5; stage++) {
        int dist = 1 << stage;
        int rot = rotations[stage];
        for (int i = 0; i < 32; i++) {
            if ((i & dist) == 0) {
                int idx1 = i;
                int idx2 = i ^ dist;
                uint64_t a = state[idx1];
                uint64_t b = state[idx2];
                a ^= b;
                b ^= rotl64(a, rot);
                state[idx1] = a;
                state[idx2] = b;
            }
        }
    }
}

void pressure_arx_scalar(uint64_t state[32]) {
    for (int c = 0; c < 8; c++) {
        uint64_t a = state[4*c];
        uint64_t b = state[4*c+1];
        uint64_t cc = state[4*c+2];
        uint64_t d = state[4*c+3];

        a += (cc ^ (cc >> 17));
        b += (d ^ (d >> 17));
        cc += (a ^ (a << 31));
        d += (b ^ (b << 31));

        state[4*c]   = a;
        state[4*c+1] = rotl64(b, 7);
        state[4*c+2] = cc;
        state[4*c+3] = rotl64(d, 19);
    }
}

void beta_iota_scalar(uint64_t state[32], int round) {
    for (int i = 0; i < 32; i++)
        state[i] ^= rc[round][i];
}

void ink_cloud_shuffle(uint64_t state[32]) {
    uint64_t temp[32];
    for (int i = 0; i < 32; i++) {
        temp[(i * 7) & 31] = rotl64(state[i], 11);
    }
    memcpy(state, temp, 32 * sizeof(uint64_t));
}

void krakken_permute_scalar_rounds(uint64_t state[32], int rounds) {
    init_rc_vectors();

    if (!state) {
        fprintf(stderr, "FATAL: krakken_permute_scalar_rounds() called with NULL state\n");
        abort();
    }
    if (rounds <= 0) return;
    if (rounds > KRAKKEN_ROUNDS) rounds = KRAKKEN_ROUNDS;

    for (int ir = 0; ir < rounds; ir++) {
        theta_scalar(state);
        tentacle_mds_scalar(state);
        rho_scalar(state);
        pi_scalar(state);
        chi_scalar(state);
        butterfly_diffusion_scalar(state);
        pressure_arx_scalar(state);
        beta_iota_scalar(state, ir);
        ink_cloud_shuffle(state);
    }
}

void krakken_permute_scalar(uint64_t state[32]) {
    krakken_permute_scalar_rounds(state, KRAKKEN_ROUNDS);
}

void krakken_hash_scalar(uint8_t *out, size_t outlen,
                         const uint8_t *in, size_t inlen) {
    if (outlen != 0 && out == NULL) {
        fprintf(stderr, "FATAL: krakken_hash_scalar() called with NULL out\n");
        abort();
    }
    if (inlen != 0 && in == NULL) {
        fprintf(stderr, "FATAL: krakken_hash_scalar() called with NULL in\n");
        abort();
    }
    if (outlen == 0) return;

    union {
        uint64_t w[32];
        uint8_t  b[256];
    } state __attribute__((aligned(32)));
    memset(state.b, 0, 256);

    const uint8_t *msg = in;
    size_t rem = inlen;

    while (rem >= 160) {
        for (int i = 0; i < 160; i++) state.b[i] ^= msg[i];
        krakken_permute_scalar(state.w);
        msg += 160; rem -= 160;
    }

    uint8_t block[160] __attribute__((aligned(32)));
    memset(block, 0, 160);
    if (rem > 0) memcpy(block, msg, rem);

    uint8_t mask = (uint8_t)(-(rem < 159));
    block[rem]  = (uint8_t)((block[rem] & ~mask) | (0x06 & mask));
    block[159]  = (uint8_t)((mask & 0x80) | ((~mask) & 0x86));
    for (int i = 0; i < 160; i++) state.b[i] ^= block[i];
    krakken_permute_scalar(state.w);

    while (outlen > 0) {
        size_t take = outlen < 160 ? outlen : 160;
        memcpy(out, state.b, take);
        out += take; outlen -= take;
        if (outlen > 0) krakken_permute_scalar(state.w);
    }

    krakken_secure_zero(state.b, 256);
    krakken_secure_zero(block,   160);
}

#ifdef KRAKKEN_MAIN
#include <sys/time.h>

static double get_time(void) {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec * 1e-6;
}

static void benchmark(void) {
    const size_t data_size = 1 * 1024 * 1024;
    const size_t hash_size = 32;
    const int    iterations = 10;
    uint8_t *data = malloc(data_size);
    uint8_t *hash = malloc(hash_size);
    if (!data || !hash) { fprintf(stderr, "malloc failed\n"); return; }

    for (size_t i = 0; i < data_size; i++)
        data[i] = (uint8_t)(i ^ (i >> 8) ^ (i >> 16) ^ (i >> 24));

    printf("Krakken-2048 Scalar Benchmark (Abyssal, %d rounds)\n", KRAKKEN_ROUNDS);
    printf("==================================================\n");
    printf("Data: %zu MB   Output: %zu bytes   Iterations: %d\n\n",
           data_size >> 20, hash_size, iterations);

    krakken_hash_scalar(hash, hash_size, data, data_size);

    double total = 0.0;
    for (int i = 0; i < iterations; i++) {
        double t0 = get_time();
        krakken_hash_scalar(hash, hash_size, data, data_size);
        total += get_time() - t0;
    }

    double avg  = total / iterations;
    double mbps = (data_size / (1024.0 * 1024.0)) / avg;
    printf("Average time : %.3f s\n", avg);
    printf("Throughput   : %.2f MB/s  (%.3f GB/s)\n", mbps, mbps / 1024.0);

    uint8_t empty[32];
    krakken_hash_scalar(empty, 32, (const uint8_t *)"", 0);
    printf("Empty hash   : ");
    for (int i = 0; i < 32; i++) printf("%02x", empty[i]);
    printf("\n");

    free(data); free(hash);
}

int main(void) {
    benchmark();
    return 0;
}
#endif