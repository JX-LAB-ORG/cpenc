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

/* MTF (move-to-front) on a BWT stream, in place: bzip2-style.  The BWT's
 * local character clustering turns into long zero runs of ranks. */
void mtf_enc_inplace(uint8_t *p, size_t n);
void mtf_dec_inplace(uint8_t *p, size_t n);

#endif /* BWT_H */
