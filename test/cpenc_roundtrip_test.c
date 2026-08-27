/* cpenc_roundtrip_test — compression stability: compress→decompress must
 * always reproduce the input exactly.  Covers sizes from 1B to 16MB,
 * special data patterns, and thousands of random payloads. */
#include "../cpenc.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int failed = 0;
#define CHK(c,m) do{if(!(c)){failed++;fprintf(stderr,"  FAIL: %s\n",m);}else{fprintf(stderr,"  OK %s\n",m);}}while(0)

static uint32_t rng_state = 0xC0FFEE;
static uint32_t rnd(void) {
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 17;
    rng_state ^= rng_state << 5;
    return rng_state;
}

/* one roundtrip: returns 1 on exact match */
static int rt(const uint8_t *in, size_t n) {
    size_t mc = cpenc_max_compressed_size(n);
    uint8_t *c = malloc(mc), *out = malloc(n);
    if (!c || !out) { free(c); free(out); return 0; }
    size_t cl = mc;
    int ok = 1;
    if (cpenc_compress(in, n, c, &cl, 2) != 0) ok = 0;
    else {
        size_t dl = n;
        if (cpenc_decompress(c, cl, out, &dl) != 0 || dl != n || memcmp(out, in, n)) ok = 0;
    }
    free(c); free(out);
    return ok;
}

int main(void) {
    printf("[c1] boundary sizes roundtrip\n");
    size_t sizes[] = {1, 2, 3, 7, 15, 16, 31, 32, 63, 64, 65, 127, 128, 129,
                      255, 256, 257, 511, 512, 1023, 1024, 4095, 4096, 4097,
                      65535, 65536, 65537, 1048576, 1048577, 16777216};
    int ok = 1;
    for (size_t si = 0; si < sizeof(sizes)/sizeof(sizes[0]); si++) {
        size_t n = sizes[si];
        uint8_t *in = malloc(n);
        for (size_t i = 0; i < n; i++) in[i] = (uint8_t)(rnd() & 0xFF);
        if (!rt(in, n)) { ok = 0; printf("  FAIL size %zu\n", n); }
        free(in);
    }
    CHK(ok, "30 boundary sizes (1B..16MB)");

    printf("[c2] special patterns roundtrip\n");
    {
        uint8_t *in = malloc(65536);
        ok = 1;
        /* all zero / all 0xFF / ramp / periodic / sparse / float-like / int4 */
        for (int p = 0; p < 8 && ok; p++) {
            for (size_t i = 0; i < 65536; i++) {
                switch (p) {
                case 0: in[i] = 0; break;
                case 1: in[i] = 0xFF; break;
                case 2: in[i] = (uint8_t)i; break;
                case 3: in[i] = (uint8_t)(i & 3); break;
                case 4: in[i] = (i % 4096) ? 0 : 0x5A; break;
                case 5: in[i] = (uint8_t)(i / 256); break;          /* float high byte */
                case 6: in[i] = (uint8_t)(i & 0x0F); break;         /* int4 */
                case 7: in[i] = (uint8_t)((i * 131) & 0xFF); break; /* periodic */
                }
            }
            if (!rt(in, 65536)) { ok = 0; printf("  FAIL pattern %d\n", p); }
        }
        CHK(ok, "8 special patterns (64KB each)");
        free(in);
    }

    printf("[c3] 2000 random roundtrips (sizes 1..128KB)\n");
    {
        ok = 1;
        for (int i = 0; i < 2000 && ok; i++) {
            size_t n = (size_t)(rnd() % 131072) + 1;
            uint8_t *in = malloc(n);
            for (size_t j = 0; j < n; j++) in[j] = (uint8_t)(rnd() & 0xFF);
            if (!rt(in, n)) { ok = 0; printf("  FAIL random #%d size %zu\n", i, n); }
            free(in);
        }
        CHK(ok, "2000 random payloads");
    }

    printf("[c4] recompression stability (compress→decompress→compress)\n");
    {
        ok = 1;
        for (int i = 0; i < 200 && ok; i++) {
            size_t n = (size_t)(rnd() % 65536) + 1;
            uint8_t *in = malloc(n);
            for (size_t j = 0; j < n; j++) in[j] = (uint8_t)(rnd() & 0xFF);
            if (!rt(in, n)) { ok = 0; }
            /* compress twice, streams must be identical (deterministic) */
            size_t mc = cpenc_max_compressed_size(n);
            uint8_t *c1 = malloc(mc), *c2 = malloc(mc);
            size_t l1 = mc, l2 = mc;
            if (cpenc_compress(in, n, c1, &l1, 2) || cpenc_compress(in, n, c2, &l2, 2) ||
                l1 != l2 || memcmp(c1, c2, l1)) { ok = 0; printf("  FAIL determinism #%d\n", i); }
            free(in); free(c1); free(c2);
        }
        CHK(ok, "200 deterministic recompressions");
    }

    printf("\n=== %d failures ===\n", failed);
    return failed ? 1 : 0;
}
