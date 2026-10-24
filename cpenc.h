#ifndef CPENC_H
#define CPENC_H
#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ─── CPEnc block header (12 + stride bytes + optional CRC) ──── */
#define CPENC_MAGIC   0x4350   /* 'CP' */
#define CPENC_VERSION 8        /* what each version added:
                                * v2: stride remainder bytes in stream;
                                * v3: 1 lag byte before orders (lag delta);
                                * v4: stride 0xFE = hex-ascii to binary layer;
                                * v5: stride 0xFD = internal segmentation;
                                * v6: stride 0xFC = BWT layer and 0xFB/0xFA =
                                *     reference layers (repeats);
                                * v7: stride 0xF9 = BWT with stored LF anchors
                                *     and RLE2; 0xFC still decodes;
                                * v8: a 4-byte CRC32 trailer (little-endian)
                                *     follows the stride payload.  The CRC is
                                *     over the ORIGINAL (pre-compression) bytes
                                *     this block reconstructs; compressed_size
                                *     includes the trailer.  Every nested block
                                *     carries its own.  A v8 block with a
                                *     missing or mismatched CRC is refused.
                                *     v1-v7 blocks have no trailer and decode
                                *     without verification. */
#define CPENC_STRIDE_RAW 0xFF  /* RAW passthrough stride */
#define CPENC_STRIDE_HEX 0xFE  /* hex-ascii text: 2:1 to binary + sep table */
#define CPENC_STRIDE_SEG 0xFD  /* segmented: per-segment lag+stride */
#define CPENC_STRIDE_BWT 0xFC  /* BWT: reorder then lag/stride (archive-mode) */
#define CPENC_STRIDE_BWT2 0xF9 /* v7 BWT: BWT + MTF + RLE2, with k stored LF
                                * anchors so decode can run k chains
                                * interleaved (see bwt.h).  Emitted when it
                                * beats the MTF-only 0xFC stream; 0xFC
                                * blocks still decode through the serial
                                * path. */
#define CPENC_STRIDE_REF 0xFB  /* reference layer: repeats as (distance,len) + patches */
#define CPENC_STRIDE_REF2 0xFA /* reference layer v2: the token stream split
                                 * into columnar sub-streams, one compressor
                                 * model each (source-code win) */
/* Reference layer: run it when this fraction of spread 16-byte windows
 * repeat, even if the head sample looked unmatchable. */
#define REF_DUP_MIN 0.10
#define CPENC_LAG_NONE 0
#define CPENC_SEG_SIZE  (512u << 10)  /* content-isolation granularity: 64KB segments spent ~35% of a structured table on per-segment tables */
#define CPENC_SEG_MIN   (1u << 20)    /* block must be >= 1MB to segment */
#define CPENC_LAG_MIN  16
#define CPENC_LAG_MAX  255      /* lag byte stores the lag itself */

typedef struct {
    uint16_t magic;           /* 0x4350 */
    uint8_t  version;         /* 1 */
    uint8_t  stride;          /* 1/2/4/8 or 0xFF = RAW */
    uint32_t original_size;
    uint32_t compressed_size; /* total after header, including per-group tANS
                               * tables and, for v8, the 4-byte CRC trailer */
} CPEncHeader;
#define CPENC_HEADER_SIZE 12
/* After the 12-byte header, there are 'stride' uint8_t order values.
 * A v8 block ends with a 4-byte little-endian CRC32 of the reconstructed
 * original bytes; it lies at header + compressed_size - 4. */

/* ─── Public API ──────────────────────────────────────────────── */

/* Compress src_len bytes from src into dst.
 *
 * Buffer contract: `*dst_len` is the capacity on input, and it must be at
 * least cpenc_max_compressed_size(src_len).  That bound is what the encoder
 * is written against; a smaller buffer is refused with -1 even when the
 * stream that would eventually be produced fits, because the encoder builds
 * and measures candidate layouts (stride, delta, BWT, reference, hex) before
 * it picks one, and those need the worst-case room.  "Exactly as large as
 * the last output for this input" is therefore not a valid size.
 *
 * dst_len: [in] capacity, [out] compressed size on success.
 * max_order: 0-4, smaller = faster but lower ratio. Tier1 uses 2, Tier2 uses 4.
 * Returns 0 on success, -1 if the buffer is below the contract or the result
 * does not fit, -2 if src_len exceeds what one block can describe
 * (0xFFFFFFFF). */
int cpenc_compress(const void *src, size_t src_len,
                   void *dst, size_t *dst_len,
                   int max_order);

/* Decompress src_len bytes from src into dst.
 *
 * dst must have room for the whole original block: a decode is not
 * incremental, and the size in the stream's header is what the caller has to
 * have believed when it sized the buffer.
 *
 * dst_len: [in] size of dst buffer, [out] actual decompressed size.
 * Returns 0 on success, -1 if dst is too small or the stream is corrupt or
 * not self-consistent, -2 if an internal allocation failed. */
int cpenc_decompress(const void *src, size_t src_len,
                     void *dst, size_t *dst_len);

/* Quick estimate of max compressed size needed (worst case). */
size_t cpenc_max_compressed_size(size_t src_len);

/* CRC-32 (IEEE) of a buffer — the same function the v8 trailer uses.  Exposed
 * so a caller that builds a block by hand (libcfs does, for RAW payloads it
 * already has in memory) can produce a valid v8 block instead of a header
 * that claims v8 and carries no checksum. */
uint32_t cpenc_crc32(const void *p, size_t n);

#ifdef __cplusplus
}
#endif
#endif /* CPENC_H */
