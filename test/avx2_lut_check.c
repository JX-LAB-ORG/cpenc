/* test/avx2_lut_check.c — differential check of the AVX2 LUT decode kernel
   against the scalar reference it replaces.

   avx2_lut16() is the x86 counterpart of ne_lut16() in cpenc_neon.h and has
   to be byte-identical to the scalar loop, or the decoder corrupts the
   stream.  Two things can go wrong independently, so both are swept:

     1. the kernel itself — every (nz, nbits) the callers can hand it, over
        random packed streams, compared against a scalar bit-cursor decode;
     2. the call-site loop — a batch stride of 2*nbits bytes and a 17-byte
        read window are only exercised by lengths that are not multiples of
        16, and a wrong stride would still pass (1) for a single batch.

   Build: cc -O2 -mavx2 test/avx2_lut_check.c -o test/avx2_lut_check
   (the Makefile target adds -mavx2 on x86-64).  On aarch64 this is a stub. */

#include <stdio.h>
#include <stdint.h>
#include <string.h>

#if defined(__x86_64__) && defined(__AVX2__)

#include "cpenc_avx2.h"

#define MAXLEN 4096
/* The kernel reads s[0..16]; the caller guarantees 17 bytes of lookahead.
   This much slack keeps the call-site guard from ever firing, so the sweep
   below really does run the vector path for every batch. */
#define PAD    64

static uint64_t rng = 0x243F6A8885A308D3ull;
static uint32_t rnd(void) {
    rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17;
    return (uint32_t)(rng >> 32);
}

/* The scalar reference: one symbol at a time off a byte-aligned bit cursor.
   Index >= 64 lands on t64[0], which is the contract avx2_lut16 mirrors
   (and what the callers get from `idx < nz ? lut[idx] : lut[0]` once t64 is
   padded with lut[0]). */
static void ref_lut16(const uint8_t *s, int nbits, const uint8_t *t64,
                      uint8_t *out) {
    uint64_t rbuf = 0; int rbc = 0; size_t p = 0;
    uint32_t mask = (uint32_t)((1u << nbits) - 1);
    for (int k = 0; k < 16; k++) {
        while (rbc < nbits) { rbuf |= (uint64_t)s[p++] << rbc; rbc += 8; }
        uint32_t idx = (uint32_t)rbuf & mask; rbuf >>= nbits; rbc -= nbits;
        out[k] = t64[idx < 64 ? idx : 0];
    }
}

/* Mirror of the tag-1 call site in cpenc.c: batch 16 symbols through the
   vector kernel while both the output and the 17-byte read window fit, then
   the scalar 4-at-a-time and 1-at-a-time tails over the same byte cursor.
   Returns how many symbols the vector path handled. */
static size_t decode_group(const uint8_t *in, size_t inlen, size_t gs, int nz,
                           int nbits, const uint8_t *lut, uint8_t *out,
                           int vec) {
    uint32_t mask = (uint32_t)((1u << nbits) - 1);
    uint64_t rbuf = 0; int rbc = 0;
    const uint8_t *p = in;
    size_t i = 0;
    if (vec) {
        uint8_t t64[64];
        for (int k = 0; k < nz; k++) t64[k] = lut[k];
        for (int k = nz; k < 64; k++) t64[k] = lut[0];
        for (; i + 16 <= gs && inlen - (size_t)(p - in) >= 17; i += 16) {
            _mm_storeu_si128((__m128i *)(out + i), avx2_lut16(p, nbits, t64));
            p += 2 * nbits;
        }
    }
    size_t vn = i;
    for (; i + 4 <= gs; i += 4) {
        while (rbc < nbits * 4 && (size_t)(p - in) < inlen) {
            rbuf |= (uint64_t)(*p++) << rbc; rbc += 8;
        }
        uint8_t i0 = (uint8_t)(rbuf & mask); rbuf >>= nbits; rbc -= nbits;
        uint8_t i1 = (uint8_t)(rbuf & mask); rbuf >>= nbits; rbc -= nbits;
        uint8_t i2 = (uint8_t)(rbuf & mask); rbuf >>= nbits; rbc -= nbits;
        uint8_t i3 = (uint8_t)(rbuf & mask); rbuf >>= nbits; rbc -= nbits;
        out[i]   = i0 < nz ? lut[i0] : lut[0];
        out[i+1] = i1 < nz ? lut[i1] : lut[0];
        out[i+2] = i2 < nz ? lut[i2] : lut[0];
        out[i+3] = i3 < nz ? lut[i3] : lut[0];
    }
    for (; i < gs; i++) {
        while (rbc < nbits && (size_t)(p - in) < inlen) {
            rbuf |= (uint64_t)(*p++) << rbc; rbc += 8;
        }
        uint8_t idx = (uint8_t)(rbuf & mask); rbuf >>= nbits; rbc -= nbits;
        out[i] = idx < nz ? lut[idx] : lut[0];
    }
    return vn;
}

static void dump(const char *what, size_t nz, int nbits, size_t gs,
                 size_t at, uint8_t got, uint8_t want) {
    fprintf(stderr, "MISMATCH %s: nz=%zu nbits=%d len=%zu first@%zu "
                    "got=%02x want=%02x\n", what, nz, nbits, gs, at, got, want);
}

int main(void) {
    static uint8_t in[MAXLEN + PAD], a[MAXLEN], b[MAXLEN];
    uint8_t lut[256], t64[64], s[64], got[16], want[16];
    size_t checks = 0, bad = 0;

    /* 1. The kernel on a single batch: every nz the callers allow, every
          nbits the kernel accepts (nbits > 6 only occurs in the synthetic
          case, but it is part of the contract the NEON kernel states). */
    for (int nz = 2; nz <= 64; nz++) {
        for (int nbits = 1; nbits <= 8; nbits++) {
            for (int rep = 0; rep < 8; rep++) {
                for (int k = 0; k < nz; k++) lut[k] = (uint8_t)rnd();
                for (int k = 0; k < nz; k++) t64[k] = lut[k];
                for (int k = nz; k < 64; k++) t64[k] = lut[0];
                for (int k = 0; k < 17; k++) s[k] = (uint8_t)rnd();
                ref_lut16(s, nbits, t64, want);
                _mm_storeu_si128((__m128i *)got, avx2_lut16(s, nbits, t64));
                checks++;
                if (memcmp(got, want, 16)) {
                    bad++;
                    if (bad <= 5) {
                        size_t k = 0;
                        while (k < 16 && got[k] == want[k]) k++;
                        dump("kernel", (size_t)nz, nbits, 16, k, got[k], want[k]);
                    }
                }
            }
        }
    }

    /* 2. The call-site loop: batch stride plus the scalar tail.  Lengths
          straddle every multiple of 16 up to 65 and then jump to a few
          thousand, which is where a mis-sized batch would show up. */
    static const size_t lens[] = {
        1, 2, 3, 4, 5, 6, 7, 8, 9, 15, 16, 17, 18, 19, 20, 31, 32, 33,
        47, 48, 49, 63, 64, 65, 127, 128, 129, 255, 256, 257,
        1000, 1023, 1024, 1025, 2048, 3000
    };
    for (int nz = 2; nz <= 64; nz++) {
        for (int nbits = 1; nbits <= 8; nbits++) {
            for (int k = 0; k < nz; k++) lut[k] = (uint8_t)rnd();
            for (unsigned li = 0; li < sizeof lens / sizeof lens[0]; li++) {
                size_t gs = lens[li];
                size_t nbytes = (gs * (size_t)nbits + 7) / 8;
                for (size_t k = 0; k < nbytes + PAD; k++) in[k] = (uint8_t)rnd();
                size_t vn = decode_group(in, nbytes + PAD, gs, nz, nbits, lut, b, 1);
                decode_group(in, nbytes + PAD, gs, nz, nbits, lut, a, 0);
                checks++;
                if (vn != (gs / 16) * 16) {
                    fprintf(stderr, "HARNESS: vector path skipped batches "
                            "(nz=%d nbits=%d len=%zu vn=%zu)\n",
                            nz, nbits, gs, vn);
                    bad++;
                } else if (memcmp(a, b, gs)) {
                    bad++;
                    if (bad <= 5) {
                        size_t k = 0;
                        while (k < gs && a[k] == b[k]) k++;
                        dump("loop", (size_t)nz, nbits, gs, k, b[k], a[k]);
                    }
                }
            }
        }
    }

    printf("avx2_lut_check: %zu cases, %zu mismatches\n", checks, bad);
    return bad != 0;
}

#else /* not x86-64 + AVX2 */

int main(void) {
#if !defined(__x86_64__)
    printf("not x86\n");
    return 0;
#else
    /* x86-64 without -mavx2 (the default ARCHFLAGS is x86-64-v2): the AVX2
       kernel is compiled out, so there is nothing to check. */
    printf("avx2_lut_check: rebuild with -mavx2 (or ARCHFLAGS=-march=x86-64-v3)\n");
    return 2;
#endif
}

#endif
