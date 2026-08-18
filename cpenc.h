#ifndef CPENC_H
#define CPENC_H
#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ─── CPEnc block header (12 + stride bytes) ─────────────────── */
#define CPENC_MAGIC   0x4350   /* 'CP' */
#define CPENC_VERSION 1
#define CPENC_STRIDE_RAW 0xFF  /* RAW passthrough stride */

typedef struct {
    uint16_t magic;           /* 0x4350 */
    uint8_t  version;         /* 1 */
    uint8_t  stride;          /* 1/2/4/8 or 0xFF = RAW */
    uint32_t original_size;
    uint32_t compressed_size; /* total after header, including per-group tANS tables */
} CPEncHeader;
#define CPENC_HEADER_SIZE 12
/* After the 12-byte header, there are 'stride' uint8_t order values */

/* ─── Public API ──────────────────────────────────────────────── */

/* Compress src_len bytes from src into dst.
 * dst_len: [in] size of dst buffer, [out] actual compressed size.
 * max_order: 0-4, smaller = faster but lower ratio. Tier1 uses 2, Tier2 uses 4.
 * Returns 0 on success, -1 if dst buffer too small, -2 on error. */
int cpenc_compress(const void *src, size_t src_len,
                   void *dst, size_t *dst_len,
                   int max_order);

/* Decompress src_len bytes from src into dst.
 * dst_len: [in] size of dst buffer, [out] actual decompressed size.
 * Returns 0 on success, -1 if dst too small or data corrupt. */
int cpenc_decompress(const void *src, size_t src_len,
                     void *dst, size_t *dst_len);

/* Quick estimate of max compressed size needed (worst case). */
size_t cpenc_max_compressed_size(size_t src_len);

#ifdef __cplusplus
}
#endif
#endif /* CPENC_H */
