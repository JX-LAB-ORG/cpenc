/* cpenc_bisect.c — shrink a failing block to the smallest length that still
 * fails to round-trip.
 *
 *   usage: cpenc_bisect FILE [ORDER]
 *
 * Compresses the first L bytes of FILE and decompresses them back, then
 * binary-searches the smallest L that fails.  Prints the stride and the head
 * of the payload at that length, which is what identifies the code path.
 */
#include "../cpenc.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static uint8_t *src;
static size_t   cap;
static uint8_t *comp, *dec;

/* 1 = round-trips, 0 = fails.  why: 0 ok, 1 compress, 2 rc, 3 length, 4 bytes */
static int try_len_why(size_t len, int order, uint8_t *stride_out, int *why, size_t *at) {
    size_t cl = cap, dl = len;
    *why = 0; *at = 0;
    if (cpenc_compress(src, len, comp, &cl, order) != 0) { *stride_out = 0xFF; *why = 1; return 0; }
    *stride_out = ((CPEncHeader *)comp)->stride;
    int rc = cpenc_decompress(comp, cl, dec, &dl);
    if (rc != 0) { *why = 2; return 0; }
    if (dl != len) { *why = 3; *at = dl; return 0; }
    for (size_t i = 0; i < len; i++)
        if (dec[i] != src[i]) { *why = 4; *at = i; return 0; }
    return 1;
}

static int try_len(size_t len, int order, uint8_t *stride_out) {
    int why; size_t at;
    return try_len_why(len, order, stride_out, &why, &at);
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s FILE [ORDER]\n", argv[0]); return 2; }
    int order = argc > 2 ? atoi(argv[2]) : 2;
    FILE *f = fopen(argv[1], "rb");
    if (!f) { perror(argv[1]); return 1; }
    fseek(f, 0, SEEK_END); long n = ftell(f); fseek(f, 0, SEEK_SET);
    src = malloc((size_t)n);
    if (fread(src, 1, (size_t)n, f) != (size_t)n) return 1;
    fclose(f);
    cap = cpenc_max_compressed_size((size_t)n);
    comp = malloc(cap); dec = malloc((size_t)n + 64);
    if (!comp || !dec) return 1;

    static const char *WHY[] = {"ok", "compress failed", "decompress returned error",
                                "length mismatch", "bytes differ"};
    uint8_t s; int why; size_t at;
    if (try_len_why((size_t)n, order, &s, &why, &at)) { printf("full length %ld round-trips\n", n); return 0; }
    printf("full length %ld fails: %s (stride %02x, at %zu)\n", n, WHY[why], s, at);

    /* smallest failing length */
    size_t lo = 1, hi = (size_t)n;
    if (try_len(lo, order, &s)) {
        while (lo + 1 < hi) {
            size_t mid = lo + (hi - lo) / 2;
            if (try_len(mid, order, &s)) lo = mid; else hi = mid;
        }
    } else {
        hi = lo;
    }
    try_len_why(hi, order, &s, &why, &at);
    printf("smallest failing length: %zu (stride %02x): %s at %zu; largest passing: %zu\n",
           hi, s, WHY[why], at, lo);

    /* what the failing block looks like */
    size_t cl = cap;
    if (cpenc_compress(src, hi, comp, &cl, order) == 0) {
        CPEncHeader *h = (CPEncHeader *)comp;
        printf("header: stride=%02x orig=%u comp=%u (total %zu)\n",
               h->stride, h->original_size, h->compressed_size, cl);
        printf("payload head:");
        for (int i = 0; i < 24 && (size_t)(CPENC_HEADER_SIZE + i) < cl; i++)
            printf(" %02x", comp[CPENC_HEADER_SIZE + i]);
        printf("\n");
    }
    /* write the failing prefix so it can be poked at directly */
    FILE *o = fopen("/root/bench/fail.bin", "wb");
    if (o) { fwrite(src, 1, hi, o); fclose(o); printf("wrote failing prefix to /root/bench/fail.bin\n"); }
    return 1;
}
