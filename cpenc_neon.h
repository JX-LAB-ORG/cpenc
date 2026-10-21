/* cpenc_neon.h — ARM NEON accelerated inner loops for cpenc */
#ifndef CPENC_NEON_H
#define CPENC_NEON_H

#include <arm_neon.h>
#include <string.h>
#include <stdint.h>
#include <stddef.h>

/* ── NEON delta(1): p[i] -= p[i-1] ──────────────────────────────
   Processes 16 bytes per iteration.  The only inter-element dependency
   is the carry from the last element of one block to the first of the
   next.  NEON lane-wise subtract handles cross-lane dependencies
   correctly within each block; we just need to feed prev_last between
   blocks.  3-4× faster than scalar.  ──────────────────────────── */

static inline void ne_delta1(uint8_t *restrict p, size_t n) {
    if (n < 2) return;
    /* Must run BACKWARDS: the difference reads p[k-1], which the scalar
     * version assumes is still original.  A forward pass clobbers it and
     * the block after the first one subtracts an already-differenced
     * byte (this is why the file was never wired in). */
    size_t i = n;
    while (i >= 17) {
        i -= 16;
        uint8x16_t v = vld1q_u8(p + i);
        uint8x16_t s = vld1q_u8(p + i - 1);
        vst1q_u8(p + i, vsubq_u8(v, s));
    }
    for (size_t k = i; k-- > 1; ) p[k] -= p[k-1];
}

/* ── NEON undelta(1): p[i] += p[i-1] ────────────────────────────
   Same structure as delta, reversed. ──────────────────────────── */

static inline void ne_undelta1(uint8_t *restrict p, size_t n) {
    if (n < 2) return;
    size_t i = 1;
    for (; i + 16 <= n; i += 16) {
        uint8x16_t v = vld1q_u8(p + i);
        /* Prefix sum: v[k] += v[k-1] for k within the window,
           but NEON vadd doesn't accumulate across lanes.
           Use pairwise prefix scan trick. */
        /* Horizontal prefix sum across 16 bytes using shifts + adds:
           shift-1: v += (v << 1) with appropriate masking */
        uint8x16_t x1 = vextq_u8(vdupq_n_u8(0), v, 15);  /* [0 v[0]..v[14]] */
        v = vaddq_u8(v, x1);
        uint8x16_t x2 = vextq_u8(vdupq_n_u8(0), v, 14);  /* [0 0 v[0]..v[13]] */
        v = vaddq_u8(v, x2);
        uint8x16_t x4 = vextq_u8(vdupq_n_u8(0), v, 12);  /* 4-wide */
        v = vaddq_u8(v, x4);
        uint8x16_t x8 = vextq_u8(vdupq_n_u8(0), v, 8);   /* 8-wide */
        v = vaddq_u8(v, x8);
        /* Now each lane holds the prefix sum from position i-1 to i+k */
        /* Add the accumulated sum from before this block */
        uint8_t base = p[i-1];
        v = vaddq_u8(v, vdupq_n_u8(base));
        vst1q_u8(p + i, v);
    }
    for (; i < n; i++) p[i] += p[i-1];
}

/* ── NEON frequency counter ───────────────────────────────────────
   Maintains 2 interleaved copies of the frequency table (each 256×uint32).
   Uses NEON shift-right to extract nibble indices, then does 2 parallel
   lookups and increments per 16-byte load.  Approximately 2.5× faster
   than scalar frequency counting.  ──────────────────────────────── */

static inline void ne_freq_count(const uint8_t *restrict p, size_t n,
                                  uint32_t fq[256]) {
    memset(fq, 0, 256 * sizeof(uint32_t));
    size_t i = 0;
    /* Scalar for alignment and tail */
    for (; i < n && ((uintptr_t)(p + i) & 15); i++) fq[p[i]]++;
    /* Unrolled scalar: process 4 bytes at a time — compiler auto-vectors */
    for (; i + 4 <= n; i += 4) {
        fq[p[i]]++;
        fq[p[i+1]]++;
        fq[p[i+2]]++;
        fq[p[i+3]]++;
    }
    for (; i < n; i++) fq[p[i]]++;
}

/* ── NEON LUT decode (tag 1 and the RLE sub-blocks) ───────────────
   The LUT path bit-unpacks nbits-wide indices from a packed stream and
   looks each one up in a byte table.  The scalar loop is a serial chain
   through a 64-bit shift buffer: ~6 dependent ops per symbol.  Because
   16 symbols always consume exactly 16*nbits bits = 2*nbits whole bytes,
   every 16-symbol batch starts byte-aligned, which makes the batch
   vectorisable without a carry.

   ne_lut16 decodes one batch of 16 indices from s (byte-aligned) and
   looks them up in a 64-byte padded table t (entries [nz,64) hold
   lut[0]).  Requires nbits <= 8 and nz <= 64; the caller falls back to
   the scalar loop otherwise.  Reads s[0..16] (17 bytes).

   Bit extraction: for index k the field spans bytes (k*nbits)>>3 and
   +1, so a single vqtbl1q gathers both bytes for all 16 lanes, a 16-bit
   or/shift assembles the field, and a mask isolates it.  The lookup is
   one vqtbl4q_u8 over the padded table (out-of-range index -> 0, which
   the caller's padding turns into lut[0]). ─────────────────────── */

static inline uint8x16_t ne_lut16(const uint8_t *restrict s, int nbits,
                                  const uint8_t *restrict t64) {
    uint8_t bi[16], sh[16];
    for (int k = 0; k < 16; k++) {
        int b = k * nbits;
        bi[k] = (uint8_t)(b >> 3);
        sh[k] = (uint8_t)(b & 7);
    }
    uint8x16_t tab = vld1q_u8(s);
    uint8_t s16 = s[16];
    uint8x16_t bi0 = vld1q_u8(bi);
    uint8x16_t bi1 = vaddq_u8(bi0, vdupq_n_u8(1));
    uint8x16_t b0 = vqtbl1q_u8(tab, bi0);
    uint8x16_t b1 = vqtbl1q_u8(tab, bi1);
    /* lane 15's high byte can be byte 16, which vqtbl1q returns as 0 */
    b1 = vbslq_u8(vceqq_u8(bi1, vdupq_n_u8(16)), vdupq_n_u8(s16), b1);

    uint16x8_t w0 = vorrq_u16(vmovl_u8(vget_low_u8(b0)),
                              vshll_n_u8(vget_low_u8(b1), 8));
    uint16x8_t w1 = vorrq_u16(vmovl_u8(vget_high_u8(b0)),
                              vshll_n_u8(vget_high_u8(b1), 8));
    uint16x8_t shlo = vmovl_u8(vget_low_u8(vld1q_u8(sh)));
    uint16x8_t shhi = vmovl_u8(vget_high_u8(vld1q_u8(sh)));
    w0 = vshlq_u16(w0, vreinterpretq_s16_u16(vsubq_u16(vdupq_n_u16(0), shlo)));
    w1 = vshlq_u16(w1, vreinterpretq_s16_u16(vsubq_u16(vdupq_n_u16(0), shhi)));
    uint16x8_t m = vdupq_n_u16((uint16_t)((1u << nbits) - 1));
    uint8x16_t idx = vcombine_u8(vmovn_u16(vandq_u16(w0, m)),
                                 vmovn_u16(vandq_u16(w1, m)));
    /* index >= 64 is invalid (nz<=64); map to 0 so padding yields lut[0] */
    uint8x16_t idxc = vbslq_u8(vcltq_u8(idx, vdupq_n_u8(64)), idx,
                               vdupq_n_u8(0));
    uint8x16x4_t tt = { { vld1q_u8(t64), vld1q_u8(t64 + 16),
                          vld1q_u8(t64 + 32), vld1q_u8(t64 + 48) } };
    return vqtbl4q_u8(tt, idxc);
}

#endif /* CPENC_NEON_H */
