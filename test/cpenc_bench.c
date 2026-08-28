/* cpenc_bench.c — compression ratio + throughput benchmark
 * Data classes: structured numerics, random, text, zeros, weights(float32),
 * byte patterns. Reports ratio and MB/s for compress/decompress. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include "../cpenc.h"

#define MB (1024*1024)
#define N  (16*MB)          /* block size per round */
#define ROUNDS 4            /* rounds for timing */

static uint64_t now_ns(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (uint64_t)ts.tv_sec * 1000000000ull + ts.tv_nsec;
}

static double mbps(uint64_t ns, size_t bytes) {
    return (double)bytes / (double)ns * 1e9 / (double)MB;
}

/* ── data generators ─────────────────────────────────────────── */
static void gen_structured_u32(uint8_t *p, size_t n) {
    uint32_t v = 0x12345678;
    for (size_t i = 0; i + 4 <= n; i += 4) {
        uint32_t *q = (uint32_t *)(p + i);
        *q = v;
        v += 0x00010003;    /* small deltas, compressible via delta+transpose */
    }
    for (size_t i = n & ~3u; i < n; i++) p[i] = (uint8_t)i;
}

static void gen_random(uint8_t *p, size_t n) {
    uint64_t s = 0x9E3779B97F4A7C15ull;
    for (size_t i = 0; i < n; i++) {
        s ^= s << 13; s ^= s >> 7; s ^= s << 17;
        p[i] = (uint8_t)(s >> 33);
    }
}

static void gen_text(uint8_t *p, size_t n) {
    static const char *words[] = {"the", "compress", "container", "disk",
        "partition", "delta", "transpose", "huffman", "stride", "memory"};
    uint64_t s = 0xDEADBEEFCAFEF00Dull;
    size_t pos = 0;
    while (pos < n) {
        const char *w = words[(s >> 40) % 10];
        s ^= s << 13; s ^= s >> 7; s ^= s << 17;
        size_t len = strlen(w);
        if (pos + len + 1 > n) break;
        memcpy(p + pos, w, len);
        pos += len;
        p[pos++] = ' ';
    }
    memset(p + pos, ' ', n - pos);
}

static void gen_zeros(uint8_t *p, size_t n) {
    memset(p, 0, n);
}

static void gen_f32_weights(uint8_t *p, size_t n) {
    /* float32 around a mean with small variance — model weight like */
    float base = 1.0f;
    uint64_t s = 0xA5A5A5A55A5A5A5Aull;
    for (size_t i = 0; i + 4 <= n; i += 4) {
        float *f = (float *)(p + i);
        *f = base + ((s >> 32) & 0xFFFF) * 1e-6f;
        s ^= s << 13; s ^= s >> 7; s ^= s << 17;
        base += 1e-5f;
    }
    for (size_t i = n & ~3u; i < n; i++) p[i] = 0;
}

static void gen_bytes(uint8_t *p, size_t n) {
    uint64_t s = 0x1111222233334444ull;
    for (size_t i = 0; i < n; i++) {
        s = s * 6364136223846793005ull + 1442695040888963407ull;
        p[i] = (uint8_t)(s >> 56);
    }
}

/* ── benchmark one class ─────────────────────────────────────── */
static void bench(const char *name, void (*gen)(uint8_t*, size_t), int max_order) {
    uint8_t *src = malloc(N), *dst = malloc(cpenc_max_compressed_size(N) + 64), *out = malloc(N);
    gen(src, N);

    /* ratio: compress whole block once */
    size_t clen = cpenc_max_compressed_size(N) + 64;
    if (cpenc_compress(src, N, dst, &clen, max_order) != 0) {
        printf("%-16s compress FAILED\n", name);
        return;
    }
    size_t rlen = N;
    if (cpenc_decompress(dst, clen, out, &rlen) != 0 || rlen != N ||
        memcmp(src, out, N) != 0) {
        printf("%-16s roundtrip MISMATCH\n", name);
        return;
    }

    /* timing */
    uint64_t t0 = now_ns();
    for (int r = 0; r < ROUNDS; r++) {
        size_t l = cpenc_max_compressed_size(N) + 64;
        cpenc_compress(src, N, dst, &l, max_order);
    }
    uint64_t tc = now_ns() - t0;
    t0 = now_ns();
    for (int r = 0; r < ROUNDS; r++) {
        size_t l = N;
        cpenc_decompress(dst, clen, out, &l);
    }
    uint64_t td = now_ns() - t0;

    printf("%-16s ratio %7.2fx   %7.2f%%   comp %6.1f MB/s   decomp %6.1f MB/s\n",
           name, (double)N / clen, (double)clen / N * 100,
           mbps(tc, (size_t)N * ROUNDS), mbps(td, (size_t)N * ROUNDS));

    free(src); free(dst); free(out);
}

int main(void) {
    printf("cpenc v%d benchmark, block=%dMB, %d rounds, order=2\n\n",
           CPENC_VERSION, (int)(N / MB), ROUNDS);
    bench("structured u32", gen_structured_u32, 2);
    bench("f32 weights", gen_f32_weights, 2);
    bench("bytes(LCG)", gen_bytes, 2);
    bench("random", gen_random, 2);
    bench("text words", gen_text, 2);
    bench("zeros", gen_zeros, 2);
    printf("\norder=4 (tier2/disk path):\n");
    bench("structured u32", gen_structured_u32, 4);
    bench("f32 weights", gen_f32_weights, 4);
    bench("text words", gen_text, 4);
    printf("\norder=0 (fastest):\n");
    bench("structured u32", gen_structured_u32, 0);
    return 0;
}
