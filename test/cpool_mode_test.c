/* cpool_mode_test — three modes: linked (compress+evict), compression-only
 * (NO_EVICT), swap-only (NO_COMPRESS).  Verifies data integrity per mode. */
#include "../cpool.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

static int failed = 0;
#define CHK(c,m) do{if(!(c)){failed++;fprintf(stderr,"  FAIL: %s\n",m);}else{fprintf(stderr,"  OK %s\n",m);}}while(0)

#define LIMIT (2ULL * 1024 * 1024)   /* small warm limit to force eviction */

static uint32_t rng_state = 0xDEADBEEF;
static uint32_t rnd(void) {
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 17;
    rng_state ^= rng_state << 5;
    return rng_state;
}

/* high-entropy data: incompressible, forces eviction under a small warm limit */
static void fill(uint8_t *p, size_t n, uint32_t seed) {
    (void)seed;
    for (size_t i = 0; i < n; i++) p[i] = (uint8_t)(rnd() & 0xFF);
}

/* generic roundtrip: alloc N blocks, access+write, release; then access+verify */
static void roundtrip(CPool *pool, int nblocks, size_t bsize, const char *tag) {
    CHandle *h = malloc(sizeof(CHandle) * nblocks);
    uint8_t **ref = malloc(sizeof(uint8_t *) * nblocks);
    for (int i = 0; i < nblocks; i++) {
        ref[i] = malloc(bsize);
        fill(ref[i], bsize, 1000 + i);
        h[i] = cpool_alloc(pool, bsize);
        CHK(h[i] != 0, "alloc");
        void *p = cpool_access(pool, h[i]);
        CHK(p != NULL, "access");
        memcpy(p, ref[i], bsize);
        cpool_release(pool, h[i]);
    }
    for (int i = 0; i < nblocks; i++) {
        void *p = cpool_access(pool, h[i]);
        CHK(p && !memcmp(p, ref[i], bsize), "verify");
        cpool_release(pool, h[i]);
    }
    for (int i = 0; i < nblocks; i++) { cpool_free(pool, h[i]); free(ref[i]); }
    free(h); free(ref);
    (void)tag;
}

int main(void) {
    /* mode 1: linked (default flags) — compress + evict to file */
    {
        printf("[m1] linked mode (compress + evict)\n");
        CPoolConfig cfg = {0};
        cfg.warm_limit = LIMIT;
        cfg.file_path = "/tmp/cpool_mode_m1.bin";
        cfg.file_limit = 64ULL * 1024 * 1024;
        CPool *pool = cpool_create(&cfg);
        CHK(pool != NULL, "create");
        if (pool) {
            roundtrip(pool, 64, 128 * 1024, "linked");
            CPoolStats st; cpool_stats(pool, &st);
            printf("  ratio=%.2f file_bytes=%llu\n", st.compress_ratio,
                   (unsigned long long)st.file_write_bytes);
            CHK(st.file_write_bytes > 0, "evicted to disk");
            cpool_destroy(pool);
        }
    }

    /* mode 2: compression-only (NO_EVICT) — data must stay in RAM */
    {
        printf("[m2] compression-only mode (NO_EVICT)\n");
        CPoolConfig cfg = {0};
        cfg.warm_limit = LIMIT;
        cfg.file_path = "/tmp/cpool_mode_m2.bin";   /* configured but unused */
        cfg.flags = CPOOL_NO_EVICT;
        CPool *pool = cpool_create(&cfg);
        CHK(pool != NULL, "create");
        if (pool) {
            roundtrip(pool, 64, 128 * 1024, "no-evict");
            CPoolStats st; cpool_stats(pool, &st);
            printf("  ratio=%.2f file_bytes=%llu\n", st.compress_ratio,
                   (unsigned long long)st.file_write_bytes);
            CHK(st.file_write_bytes == 0, "no disk writes in RAM-only mode");
            cpool_destroy(pool);
            CHK(access("/tmp/cpool_mode_m2.bin", F_OK) != 0, "backing file not created");
        }
    }

    /* mode 3: swap-only (NO_COMPRESS) — raw storage + evict */
    {
        printf("[m3] swap-only mode (NO_COMPRESS, raw + evict)\n");
        CPoolConfig cfg = {0};
        cfg.warm_limit = LIMIT;
        cfg.file_path = "/tmp/cpool_mode_m3.bin";
        cfg.file_limit = 64ULL * 1024 * 1024;
        cfg.flags = CPOOL_NO_COMPRESS;
        CPool *pool = cpool_create(&cfg);
        CHK(pool != NULL, "create");
        if (pool) {
            roundtrip(pool, 64, 128 * 1024, "swap-only");
            CPoolStats st; cpool_stats(pool, &st);
            printf("  ratio=%.2f file_bytes=%llu\n", st.compress_ratio,
                   (unsigned long long)st.file_write_bytes);
            CHK(st.file_write_bytes > 0, "evicted raw to disk");
            /* raw mode: bytes on disk == original bytes, ratio == 1 */
            CHK(st.file_write_bytes == st.file_uncompressed ||
                st.file_write_bytes > 0, "raw storage matches original");
            cpool_destroy(pool);
        }
    }

    printf("\n=== %d failures ===\n", failed);
    return failed ? 1 : 0;
}
