/* cpenc_avx2.h — x86-64 AVX2 accelerated inner loops for cpenc */
#ifndef CPENC_AVX2_H
#define CPENC_AVX2_H

/* Self-guarding: on anything that is not x86-64 with AVX2 enabled the
   header expands to nothing, so the callers' scalar loops stay in charge
   and the compressed bitstream is bit-identical everywhere.  <immintrin.h>
   is only ever pulled in when the target really has the intrinsics. */
#if defined(__x86_64__) && defined(__AVX2__)

#include <immintrin.h>
#include <stdint.h>
#include <stddef.h>

/* ── AVX2 LUT decode (tag 1 and the RLE sub-blocks) ───────────────
   Same input/output contract as ne_lut16() in cpenc_neon.h: decode one
   batch of 16 symbols from a byte-aligned packed stream s, look each one
   up in the 64-byte padded table t64, return the 16 result bytes.

   nbits in [1,8] and nz <= 64; the caller falls back to scalar otherwise.
   Reads s[0..16] (17 bytes, the last one only for the highest symbol's
   upper byte).  An index >= 64 decodes to t64[0] — that is what NEON's
   vqtbl4q out-of-range rule yields and what the scalar `idx < nz ? lut[idx]
   : lut[0]` yields once t64 is padded with lut[0].

   Why batching is legal at all: 16 symbols always consume exactly
   16*nbits bits = 2*nbits whole bytes, so every batch starts byte-aligned
   and needs no carry from the previous one.  Within a batch symbol k sits
   at bit k*nbits, so with nbits <= 8 its field spans at most the two bytes
   (k*nbits)>>3 and +1 — one 16-bit little-endian window per lane.

   AVX2 has no cross-lane byte shuffle and no 64-byte table instruction, so
   the two things NEON does in one op each cost more here: the byte gather
   is a vpshufb (which reads only the low nibble of its index) plus an
   explicit fix-up for the one lane whose second byte is byte 16, and the
   64-entry table read is four 16-byte vpshufb lookups selected by the
   index's high bits.  The per-lane variable shift is vpsrlvd, which only
   exists for 32-bit lanes, hence the widen/shift/narrow.  ────────────── */

/* Field byte offset (k*nbits)>>3 and in-byte shift (k*nbits)&7 for
   symbols k=0..15, indexed by nbits-1.  Precomputed because the arithmetic
   costs more vector ops than the decode it feeds. */
static const uint8_t avx2_lut_bi[8][16] = {
    { 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1 },
    { 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3 },
    { 0, 0, 0, 1, 1, 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5 },
    { 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7 },
    { 0, 0, 1, 1, 2, 3, 3, 4, 5, 5, 6, 6, 7, 8, 8, 9 },
    { 0, 0, 1, 2, 3, 3, 4, 5, 6, 6, 7, 8, 9, 9, 10, 11 },
    { 0, 0, 1, 2, 3, 4, 5, 6, 7, 7, 8, 9, 10, 11, 12, 13 },
    { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 },
};
static const uint8_t avx2_lut_sh[8][16] = {
    { 0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7 },
    { 0, 2, 4, 6, 0, 2, 4, 6, 0, 2, 4, 6, 0, 2, 4, 6 },
    { 0, 3, 6, 1, 4, 7, 2, 5, 0, 3, 6, 1, 4, 7, 2, 5 },
    { 0, 4, 0, 4, 0, 4, 0, 4, 0, 4, 0, 4, 0, 4, 0, 4 },
    { 0, 5, 2, 7, 4, 1, 6, 3, 0, 5, 2, 7, 4, 1, 6, 3 },
    { 0, 6, 4, 2, 0, 6, 4, 2, 0, 6, 4, 2, 0, 6, 4, 2 },
    { 0, 7, 6, 5, 4, 3, 2, 1, 0, 7, 6, 5, 4, 3, 2, 1 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static inline __m128i avx2_lut16(const uint8_t *restrict s, int nbits,
                                 const uint8_t *restrict t64) {
    __m128i bi = _mm_loadu_si128((const __m128i *)avx2_lut_bi[nbits - 1]);
    __m128i sh = _mm_loadu_si128((const __m128i *)avx2_lut_sh[nbits - 1]);

    /* Gather both bytes of every 16-bit field.  vpshufb indexes within its
       128-bit source and uses only the low nibble, so bi (0..15) is exact;
       bi+1 reaches 16 for the top symbol, which the shuffle folds back to
       tab[0], so that lane is patched from s[16] explicitly. */
    __m128i tab = _mm_loadu_si128((const __m128i *)s);
    __m128i bi1 = _mm_add_epi8(bi, _mm_set1_epi8(1));
    __m128i b0 = _mm_shuffle_epi8(tab, bi);
    __m128i b1 = _mm_shuffle_epi8(tab, bi1);
    b1 = _mm_blendv_epi8(b1, _mm_set1_epi8((char)s[16]),
                         _mm_cmpeq_epi8(bi1, _mm_set1_epi8(16)));

    /* Assemble the little-endian field words, then shift each down by its
       own bit offset.  vpsrlvd is 32-bit only, so widen, shift, narrow. */
    __m128i wlo = _mm_unpacklo_epi8(b0, b1);
    __m128i whi = _mm_unpackhi_epi8(b0, b1);
    __m256i m = _mm256_set1_epi32((1 << nbits) - 1);
    __m256i x0 = _mm256_and_si256(_mm256_srlv_epi32(
        _mm256_cvtepu16_epi32(wlo), _mm256_cvtepu8_epi32(sh)), m);
    __m256i x1 = _mm256_and_si256(_mm256_srlv_epi32(
        _mm256_cvtepu16_epi32(whi),
        _mm256_cvtepu8_epi32(_mm_srli_si128(sh, 8))), m);
    __m128i idx = _mm_packus_epi16(
        _mm_packus_epi32(_mm256_castsi256_si128(x0),
                         _mm256_extracti128_si256(x0, 1)),
        _mm_packus_epi32(_mm256_castsi256_si128(x1),
                         _mm256_extracti128_si256(x1, 1)));

    /* An index >= 64 is off the table (nz <= 64): fold it to 0 so the read
       lands on t64[0] = lut[0].  Also keeps bit 7 clear, which vpshufb
       would otherwise read as "return zero". */
    __m128i in64 = _mm_cmpeq_epi8(_mm_subs_epu8(idx, _mm_set1_epi8(63)),
                                  _mm_setzero_si128());
    idx = _mm_and_si128(idx, in64);

    /* 64 entries = four 16-byte sub-tables chosen by index bits 4-5. */
    __m128i sel = _mm_and_si128(_mm_srli_epi16(idx, 4), _mm_set1_epi8(3));
    __m128i lo = _mm_cmpeq_epi8(_mm_and_si128(sel, _mm_set1_epi8(1)),
                                _mm_set1_epi8(1));
    __m128i hi = _mm_cmpeq_epi8(_mm_and_si128(sel, _mm_set1_epi8(2)),
                                _mm_set1_epi8(2));
    __m128i r0 = _mm_shuffle_epi8(_mm_loadu_si128((const __m128i *)(t64 + 0)), idx);
    __m128i r1 = _mm_shuffle_epi8(_mm_loadu_si128((const __m128i *)(t64 + 16)), idx);
    __m128i r2 = _mm_shuffle_epi8(_mm_loadu_si128((const __m128i *)(t64 + 32)), idx);
    __m128i r3 = _mm_shuffle_epi8(_mm_loadu_si128((const __m128i *)(t64 + 48)), idx);
    return _mm_blendv_epi8(_mm_blendv_epi8(r0, r1, lo),
                           _mm_blendv_epi8(r2, r3, lo), hi);
}

#endif /* __x86_64__ && __AVX2__ */

#endif /* CPENC_AVX2_H */
