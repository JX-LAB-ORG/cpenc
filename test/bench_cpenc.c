/* bench_cpenc.c — codec-level benchmark over a real file.
 *
 * Reports ratio and compress/decompress throughput for one block size and
 * order, plus a histogram of which stride (inner model) cpenc picked per
 * block — that is the part of "how well does it do" that a single ratio
 * number hides: it says whether the data went through the reference layer,
 * the BWT layer, plain stride/delta, or came out as RAW.
 *
 *   usage: bench_cpenc FILE [BLOCK_KB] [ORDER] [REPS]
 *
 * The three phases are strictly separate: compress every block (timed),
 * decompress those buffers (timed), then verify (untimed).  Compressing
 * inside the decompression loop — as an earlier version of this file did —
 * reports compress+decompress as "D MB/s", and the giveaway is C and D
 * coming out nearly equal.
 */
#include "../cpenc.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MB (1024 * 1024)

static double now_s(void) {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (double)t.tv_sec + (double)t.tv_nsec / 1e9;
}

static const char *stride_name(uint8_t s) {
    switch (s) {
    case CPENC_STRIDE_RAW:  return "RAW";
    case CPENC_STRIDE_HEX:  return "HEX";
    case CPENC_STRIDE_SEG:  return "SEG";
    case CPENC_STRIDE_BWT:  return "BWT";
    case CPENC_STRIDE_REF:  return "REF";
    case CPENC_STRIDE_REF2: return "REF2";
    default:                return "STRIDE";
    }
}

static int stride_slot(uint8_t s) {
    switch (s) {
    case CPENC_STRIDE_RAW:  return 0;
    case CPENC_STRIDE_HEX:  return 1;
    case CPENC_STRIDE_SEG:  return 2;
    case CPENC_STRIDE_BWT:  return 3;
    case CPENC_STRIDE_REF:  return 4;
    case CPENC_STRIDE_REF2: return 5;
    default:                return 0;
    }
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "usage: %s FILE [BLOCK_KB] [ORDER] [REPS]\n", argv[0]);
        return 2;
    }
    const char *path = argv[1];
    size_t blk = (argc > 2 ? (size_t)atoi(argv[2]) : 8192) * 1024;
    int order  = (argc > 3 ? atoi(argv[3]) : 2);
    int reps   = (argc > 4 ? atoi(argv[4]) : 3);
    if (blk == 0 || reps < 1) { fprintf(stderr, "bad args\n"); return 2; }

    FILE *f = fopen(path, "rb");
    if (!f) { perror(path); return 1; }
    fseek(f, 0, SEEK_END);
    long n = ftell(f);
    fseek(f, 0, SEEK_SET);
    if (n <= 0) { fprintf(stderr, "empty file\n"); return 1; }
    uint8_t *src = malloc((size_t)n);
    if (!src || fread(src, 1, (size_t)n, f) != (size_t)n) { fprintf(stderr, "read\n"); return 1; }
    fclose(f);

    size_t nblk = ((size_t)n + blk - 1) / blk;
    uint8_t **cbuf = calloc(nblk, sizeof(*cbuf));
    size_t   *clen = calloc(nblk, sizeof(*clen));
    size_t   *olen = calloc(nblk, sizeof(*olen));
    uint8_t  *dec  = malloc(blk);
    if (!cbuf || !clen || !olen || !dec) { fprintf(stderr, "alloc\n"); return 1; }

    size_t bound = cpenc_max_compressed_size(blk);
    /* One scratch buffer, and each block keeps only its actual compressed
     * size.  Allocating a full worst-case buffer per block — as this did —
     * costs about 1.5x the input, which is fine for a 64MB corpus and fatal
     * for the 673MB disk image this has to measure. */
    uint8_t *scratch = malloc(bound);
    if (!scratch) { fprintf(stderr, "alloc scratch\n"); return 1; }
    for (size_t i = 0; i < nblk; i++) cbuf[i] = NULL;   /* sized after compressing */

    uint64_t hist[6] = {0};
    size_t total_orig = 0, total_comp = 0;
    double tc_best = 1e18, td_best = 1e18;

    for (int r = 0; r < reps; r++) {
        /* ── phase 1: compress every block ── */
        memset(hist, 0, sizeof(hist));
        size_t co = 0, cc = 0;
        double t0 = now_s();
        for (size_t i = 0; i < nblk; i++) {
            size_t off = i * blk;
            size_t len = (size_t)n - off < blk ? (size_t)n - off : blk;
            size_t cl = bound;
            if (cpenc_compress(src + off, len, scratch, &cl, order) != 0) {
                fprintf(stderr, "compress fail at block %zu\n", i);
                return 1;
            }
            /* keep only this block's actual bytes: the timing loop runs
             * more than once, so the previous copy goes first */
            free(cbuf[i]);
            cbuf[i] = malloc(cl ? cl : 1);
            if (!cbuf[i]) { fprintf(stderr, "alloc block %zu\n", i); return 1; }
            memcpy(cbuf[i], scratch, cl);
            clen[i] = cl;
            olen[i] = len;
            co += len; cc += cl;
            hist[stride_slot(((CPEncHeader *)cbuf[i])->stride)]++;
        }
        double tc = now_s() - t0;

        /* ── phase 2: decompress those buffers ── */
        double t1 = now_s();
        for (size_t i = 0; i < nblk; i++) {
            size_t dl = blk;
            if (cpenc_decompress(cbuf[i], clen[i], dec, &dl) != 0) {
                fprintf(stderr, "decompress fail at block %zu\n", i);
                return 1;
            }
        }
        double td = now_s() - t1;

        if (tc < tc_best) tc_best = tc;
        if (td < td_best) td_best = td;
        total_orig = co; total_comp = cc;
    }

    /* ── phase 3: verify (untimed) ── */
    int bad = 0;
    for (size_t i = 0; i < nblk; i++) {
        size_t off = i * blk;
        size_t dl = blk;
        if (cpenc_decompress(cbuf[i], clen[i], dec, &dl) != 0
            || dl != olen[i] || memcmp(dec, src + off, olen[i]) != 0) bad++;
    }

    uint64_t nb = 0;
    for (int i = 0; i < 6; i++) nb += hist[i];

    const char *base = strrchr(path, '/');
    base = base ? base + 1 : path;
    printf("%-12s blk=%-5zuK ord=%d  %7.2f -> %7.2f MB  ratio %6.2fx  "
           "C %7.2f MB/s  D %8.2f MB/s  [",
           base, blk / 1024, order,
           (double)total_orig / MB, (double)total_comp / MB,
           (double)total_orig / (double)(total_comp ? total_comp : 1),
           (double)total_orig / MB / tc_best, (double)total_orig / MB / td_best);
    int first = 1;
    for (int i = 0; i < 6; i++) {
        if (!hist[i]) continue;
        const char *nm = stride_name(i == 0 ? CPENC_STRIDE_RAW : i == 1 ? CPENC_STRIDE_HEX
                                 : i == 2 ? CPENC_STRIDE_SEG : i == 3 ? CPENC_STRIDE_BWT
                                 : i == 4 ? CPENC_STRIDE_REF : CPENC_STRIDE_REF2);
        printf("%s%s %.0f%%", first ? "" : " ", nm, 100.0 * (double)hist[i] / (double)nb);
        first = 0;
    }
    printf("]%s\n", bad ? "  ROUNDTRIP-FAIL" : "");

    for (size_t i = 0; i < nblk; i++) free(cbuf[i]);
    free(cbuf); free(clen); free(olen); free(dec); free(scratch); free(src);
    return bad ? 1 : 0;
}
