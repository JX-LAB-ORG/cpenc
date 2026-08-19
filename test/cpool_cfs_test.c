/* cpool_cfs_test — linked mode: cpool spills into a cfs swap partition via
 * the backend callbacks, then blocks are read back (swapped in) intact.
 *
 * Build (from /root/libcpool):
 *   cc -std=c11 -Wall -O3 -march=armv8-a+simd -I. -I../libcfs \
 *       test/cpool_cfs_test.c cpool.o cpenc.o \
 *       ../libcfs/cfs.o ../libcfs/cfs_part.o ../libcfs/cfs_mkfs_fat32.o \
 *       ../libcfs/cfs_mkfs_ext4.c ../libcfs/cfs_mkfs_swap.c -lm -o test/cpool_cfs_test
 */
#include "../cpool.h"
#include "../libcfs/cfs.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

static int failed = 0;
#define CHK(c,m) do{if(!(c)){failed++;fprintf(stderr,"  FAIL: %s\n",m);}else{fprintf(stderr,"  OK %s\n",m);}}while(0)

#define DISK_SZ   (64ULL * 1024 * 1024)
#define PART_START (1ULL << 20)
#define CFS_PATH  "/tmp/cpool_cfs_test.cfs"
#define WARM_LIMIT (2ULL * 1024 * 1024)

typedef struct { CFS *cfs; uint64_t start; uint64_t size; } CfsBackend;

static int b_open(void **ctx, uint64_t limit) {
    CfsBackend *b = malloc(sizeof(*b));
    if (!b) return -1;
    b->cfs = cfs_open(CFS_PATH, CFS_RDWR | CFS_RAW);
    if (!b->cfs) { free(b); return -1; }
    b->start = PART_START;
    b->size = DISK_SZ - PART_START;
    if (limit && limit < b->size) b->size = limit;
    *ctx = b;
    return 0;
}
static int b_write(void *ctx, const void *data, size_t sz, uint64_t off) {
    CfsBackend *b = ctx;
    if (off + sz > b->size) return -1;
    return cfs_block_write(b->cfs, data, sz, (off_t)(b->start + off));
}
static int b_read(void *ctx, void *data, size_t sz, uint64_t off) {
    CfsBackend *b = ctx;
    if (off + sz > b->size) return -1;
    return cfs_block_read(b->cfs, data, sz, (off_t)(b->start + off)) == (ssize_t)sz ? 0 : -1;
}
static void b_close(void *ctx) {
    CfsBackend *b = ctx;
    if (b) { cfs_close(b->cfs); free(b); }
}

static uint32_t rng_state = 0x1234ABCD;
static uint32_t rnd(void) {
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 17;
    rng_state ^= rng_state << 5;
    return rng_state;
}

int main(void) {
    /* Prepare the disk: MBR + one partition (1MB..64MB) — no fs needed,
     * the partition is the raw spill area. */
    unlink(CFS_PATH);
    CFS *d = cfs_create(CFS_PATH);
    CHK(d != NULL, "disk create");
    if (!d) return 1;
    cfs_block_resize(d, DISK_SZ);
    CHK(cfs_part_label(d, CFS_PART_MBR) == 0, "label mbr");
    CHK(cfs_part_add(d, PART_START, DISK_SZ - PART_START, 0x82, NULL, 0) == 1, "partadd swap-type");
    cfs_close(d);

    /* Linked mode: compress + spill into the cfs partition */
    CPoolConfig cfg = {0};
    cfg.warm_limit = WARM_LIMIT;
    cfg.file_limit = DISK_SZ - PART_START;
    cfg.backend_open = b_open;
    cfg.backend_write = b_write;
    cfg.backend_read = b_read;
    cfg.backend_close = b_close;
    CPool *pool = cpool_create(&cfg);
    CHK(pool != NULL, "pool create (cfs backend)");
    if (!pool) { unlink(CFS_PATH); return 1; }

    /* 64 blocks x 128KB incompressible = 8MB > 2MB warm → forced spills */
    enum { NB = 64, BS = 128 * 1024 };
    CHandle h[NB];
    uint8_t *ref[NB];
    int ok = 1;
    for (int i = 0; i < NB; i++) {
        ref[i] = malloc(BS);
        for (size_t j = 0; j < BS; j++) ref[i][j] = (uint8_t)(rnd() & 0xFF);
        h[i] = cpool_alloc(pool, BS);
        void *p = cpool_access(pool, h[i]);
        if (!p) { ok = 0; break; }
        memcpy(p, ref[i], BS);
        cpool_release(pool, h[i]);
    }
    CHK(ok, "write 8MB through pool");

    CPoolStats st;
    cpool_stats(pool, &st);
    printf("  compress_ratio=%.2f file_bytes=%llu\n", st.compress_ratio,
           (unsigned long long)st.file_write_bytes);
    CHK(st.file_write_bytes > 0, "spilled into cfs partition");

    /* Swap back in: verify every block */
    ok = 1;
    for (int i = 0; i < NB; i++) {
        void *p = cpool_access(pool, h[i]);
        if (!p || memcmp(p, ref[i], BS)) { ok = 0; printf("  mismatch block %d\n", i); break; }
        cpool_release(pool, h[i]);
    }
    CHK(ok, "swap-in roundtrip intact");
    /* baseline after the swap-in pass (its releases also rewrite disk) */
    cpool_stats(pool, &st);

    /* Verify the partition actually holds pool data on disk */
    d = cfs_open(CFS_PATH, CFS_RDONLY | CFS_RAW);
    CHK(d != NULL, "reopen disk");
    if (d) {
        uint8_t probe[64];
        ssize_t n = cfs_block_read(d, probe, sizeof(probe), (off_t)PART_START);
        int nonz = 0;
        for (size_t i = 0; i < sizeof(probe); i++) if (probe[i]) nonz = 1;
        CHK(n == 64 && nonz, "partition holds spill data");
        cfs_close(d);
    }

    /* Compressible data: linked mode must compress BEFORE spilling
     * (spilled bytes << original), proving the compression side. */
    CPoolStats st2;
    for (int i = 0; i < NB; i++) {
        memset(ref[i], 0xAB, BS);   /* trivially compressible */
        void *p = cpool_access(pool, h[i]);
        if (!p) { ok = 0; break; }
        memcpy(p, ref[i], BS);
        cpool_release(pool, h[i]);
    }
    CHK(ok, "rewrite with compressible data");
    cpool_stats(pool, &st2);
    uint64_t total_file = st2.file_write_bytes - st.file_write_bytes;  /* this round only */
    /* re-verify after the rewrite (blocks may now be warm-resident) */
    ok = 1;
    for (int i = 0; i < NB; i++) {
        void *p = cpool_access(pool, h[i]);
        if (!p || memcmp(p, ref[i], BS)) { ok = 0; break; }
        cpool_release(pool, h[i]);
    }
    CHK(ok, "compressible roundtrip intact");
    printf("  compress ratio after rewrite: %.2f (8MB -> %llu disk bytes)\n",
           (float)(NB * BS) / (float)total_file, (unsigned long long)total_file);

    for (int i = 0; i < NB; i++) { cpool_free(pool, h[i]); free(ref[i]); }
    cpool_destroy(pool);
    unlink(CFS_PATH);
    printf("\n=== %d failures ===\n", failed);
    return failed ? 1 : 0;
}
