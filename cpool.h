#ifndef CPOOL_H
#define CPOOL_H
#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ─── Opaque handle ────────────────────────────────────────── */
typedef uint64_t CHandle;

/* ─── Configuration ────────────────────────────────────────── */
/* Mode flags: 0 (default) = compress + evict (linked mode).
 *   CPOOL_NO_COMPRESS : store blocks raw (swap-only mode)
 *   CPOOL_NO_EVICT    : RAM-only (compression-only mode, no disk spill) */
#define CPOOL_NO_COMPRESS (1u << 0)
#define CPOOL_NO_EVICT    (1u << 1)

typedef struct {
    size_t warm_limit;            /* max compressed RAM bytes */
    const char *file_path;        /* disk backing file (NULL=RAM only) */
    size_t file_limit;            /* max disk bytes, 0=unlimited */
    int    compress_cores[8];     /* CPU affinity for compress workers */
    int    n_compress_cores;      /* 0=no pinning */
    int    batch_flush_size;      /* bytes to batch before disk flush */
    unsigned flags;               /* CPOOL_NO_* bits; 0 = linked mode */

    /* Optional custom backing backend (e.g. a cfs partition).
     * When backend_open is set, file_path is ignored. */
    int  (*backend_open)(void **ctx, uint64_t limit);
    int  (*backend_write)(void *ctx, const void *data, size_t sz, uint64_t off);
    int  (*backend_read)(void *ctx, void *data, size_t sz, uint64_t off);
    void (*backend_close)(void *ctx);
} CPoolConfig;

/* ─── Opaque pool handle ───────────────────────────────────── */
typedef struct CPool CPool;

/* ─── Lifecycle ────────────────────────────────────────────── */
CPool *cpool_create(const CPoolConfig *cfg);
void   cpool_destroy(CPool *pool);

/* ─── Alloc / Free ─────────────────────────────────────────── */
CHandle cpool_alloc(CPool *pool, size_t size);
void    cpool_free(CPool *pool, CHandle h);

/* ─── Access / Release (read / write temporary buffer) ─────── */
void  *cpool_access(CPool *pool, CHandle h);
void   cpool_release(CPool *pool, CHandle h);

/* ─── Hint (optional) ──────────────────────────────────────── */
void   cpool_hint_cold(CPool *pool, CHandle h);

/* ─── Statistics ───────────────────────────────────────────── */
typedef struct {
    size_t   warm_used;           /* compressed bytes in RAM */
    size_t   warm_uncompressed;   /* equivalent uncompressed */
    size_t   file_used;           /* bytes on disk */
    size_t   file_uncompressed;
    float    compress_ratio;      /* overall */
    uint64_t compress_count;
    uint64_t decompress_count;
    uint64_t file_write_bytes;
    uint64_t file_read_bytes;
} CPoolStats;
void cpool_stats(CPool *pool, CPoolStats *out);

/* ─── Internals (for test/debug) ───────────────────────────── */
size_t cpool_warm_free(const CPool *pool);

#ifdef __cplusplus
}
#endif
#endif /* CPOOL_H */
