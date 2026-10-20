/* cpenc_rle_test.c — regression for the tag-4 (RLE) run-stream bound.
 *
 * The RLE tag stores runs as a byte stream: a run of 255 or more is written
 * as the byte 255 plus a LEB128 remainder, so the stream is LONGER than the
 * run count.  The decoder used to size its run buffer (and reject the block)
 * on gs/2+2, which is the bound on the run COUNT given the encoder's
 * zero-density gate — not on the encoded byte count.  A block whose run
 * stream came out a few bytes longer was rejected as corrupt, which is what
 * made a 256KB BWT block of libLLVM fail to decompress.
 *
 * The shapes below reproduce that structure: zero blocks and non-zero blocks
 * alternating, with the zero blocks long enough to need the 255+LEB128 escape
 * and the non-zero blocks byte-by-byte (so the run count sits at its
 * maximum, nv+1).  Every case must round-trip.
 */
#include "../cpenc.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int failed = 0;
#define CHK(c, m) do { if (!(c)) { failed++; fprintf(stderr, "  FAIL: %s\n", m); } \
                       else fprintf(stderr, "  OK %s\n", m); } while (0)

static uint64_t rng = 0x9E3779B97F4A7C15ull;
static uint8_t nz(void) {           /* a non-zero byte */
    rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17;
    return (uint8_t)((rng >> 33) | 1);
}

/* nseg segments of [zlen zeros][nlen non-zero bytes] */
static size_t build(uint8_t *p, int nseg, size_t zlen, size_t nlen) {
    size_t o = 0;
    for (int s = 0; s < nseg; s++) {
        memset(p + o, 0, zlen); o += zlen;
        for (size_t k = 0; k < nlen; k++) p[o++] = nz();
    }
    return o;
}

static void roundtrip(const char *what, const uint8_t *src, size_t n, int order) {
    size_t cap = cpenc_max_compressed_size(n);
    uint8_t *comp = malloc(cap), *dec = malloc(n);
    size_t cl = cap, dl = n;
    if (cpenc_compress(src, n, comp, &cl, order) != 0) { CHK(0, what); free(comp); free(dec); return; }
    int rc = cpenc_decompress(comp, cl, dec, &dl);
    uint8_t stride = ((CPEncHeader *)comp)->stride;
    CHK(rc == 0 && dl == n && memcmp(dec, src, n) == 0, what);
    if (rc != 0 || dl != n || memcmp(dec, src, n) != 0)
        fprintf(stderr, "       (n=%zu stride=%02x cl=%zu rc=%d dl=%zu)\n", n, stride, cl, rc, dl);
    free(comp); free(dec);
}

int main(void) {
    /* The failing case had zeros just over half the group, ~131k runs, and
     * 20 runs long enough to escape.  zlen = nlen+1 keeps the gate satisfied
     * (zeros > gs/2) while leaving the run count at its maximum. */
    static uint8_t buf[8 << 20];
    struct { int nseg; size_t zlen, nlen; } cases[] = {
        { 20, 6501, 6500 },     /* 260020 bytes: the shape that failed */
        { 20, 255,  254  },     /* smallest escape */
        { 4,  32768, 32767 },   /* one long escape per segment */
        { 40, 1023, 1022 },
        { 2,  131071, 131070 },
    };
    for (unsigned i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
        size_t n = build(buf, cases[i].nseg, cases[i].zlen, cases[i].nlen);
        char m[128];
        for (int o = 0; o <= 2; o += 2) {
            snprintf(m, sizeof(m), "RLE bound: %d x [%zu zeros,%zu nz] = %zu bytes, order %d",
                     cases[i].nseg, cases[i].zlen, cases[i].nlen, n, o);
            roundtrip(m, buf, n, o);
        }
    }
    printf("\n=== %d failures ===\n", failed);
    return failed ? 1 : 0;
}
