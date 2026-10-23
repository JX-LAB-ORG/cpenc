#ifndef BWT_H
#define BWT_H
#include <stddef.h>
#include <stdint.h>

/* Burrows-Wheeler transform (self-contained, doubling suffix sort with
 * 2-pass radix on rank pairs; decode via LF mapping, O(n)).
 *
 * bwt_enc: output[n] = BWT of s[n], returns the index of the original
 *          string's first row (needed for decode).
 * bwt_dec: bwt_dec(bwt[n], n, orig_idx, out[n]) restores s exactly.
 *          Returns 0, or -1 when the stream is not a valid BWT: the primary
 *          index is out of range or the LF chain leaves the block.  The
 *          decoder feeds this function bytes from an untrusted container, and
 *          following a corrupt chain reads outside both buffers. */
size_t bwt_enc(const uint8_t *s, size_t n, uint8_t *out);
int    bwt_dec(const uint8_t *bwt, size_t n, size_t orig_idx, uint8_t *out);

/* v7 anchored decode: the LF cycle is one serial dependency chain, and its
 * latency is what limits BWT decode.  Storing k positions along the cycle
 * lets k chains advance concurrently (single thread, memory-level
 * parallelism).  n1 is the row count including the EOF row.
 *
 * bwt_make_anchors: after bwt_enc, record anchor[j] = the LF position after
 *   j*m steps, m = n1/k.  anchors must hold k entries; k >= 1.
 * bwt_dec_anchors: rebuild the same LF as bwt_dec and run the k chains
 *   interleaved.  returns 0, or -1 for an invalid stream.  The anchors come
 *   from an untrusted container, so they are validated: 8 <= k <= 16,
 *   n1 >= k, anchors[0] == 0 (the cycle start), every anchor < n1, and all
 *   k pairwise distinct.  Every LF value is range-checked as the table is
 *   built, so no position a chain can follow lies outside the block;
 *   between that and the anchor checks, no step can read outside bwt[] or
 *   write outside out[]. */
int bwt_make_anchors(const uint8_t *bwt, size_t n1, size_t orig_idx,
                     uint32_t *anchors, int k);
int bwt_dec_anchors(const uint8_t *bwt, size_t n1, size_t orig_idx,
                    const uint32_t *anchors, int k, uint8_t *out);

/* MTF (move-to-front) on a BWT stream, in place: bzip2-style.  The BWT's
 * local character clustering turns into long zero runs of ranks. */
void mtf_enc_inplace(uint8_t *p, size_t n);
void mtf_dec_inplace(uint8_t *p, size_t n);

#endif /* BWT_H */
