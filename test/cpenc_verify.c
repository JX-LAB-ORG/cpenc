/* cpenc_verify.c — block-size sweep that checks round-trip, and on failure
 * reports enough to reproduce: block size, offset, length, and the stride the
 * codec picked (which says which inner model to look at).
 *
 *   usage: cpenc_verify FILE [ORDER] [BLK_KB...]
 */
#include "../cpenc.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const char *sname(uint8_t s) {
    switch (s) {
    case CPENC_STRIDE_RAW:  return "RAW";
    case CPENC_STRIDE_HEX:  return "HEX";
    case CPENC_STRIDE_SEG:  return "SEG";
    case CPENC_STRIDE_BWT:  return "BWT";
    case CPENC_STRIDE_REF:  return "REF";
    case CPENC_STRIDE_REF2: return "REF2";
    default:                return "?";
    }
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s FILE [ORDER] [BLK_KB...]\n", argv[0]); return 2; }
    int order = argc > 2 ? atoi(argv[2]) : 2;
    size_t blks[16];
    int nblk = 0;
    if (argc > 3) {
        for (int i = 3; i < argc && nblk < 16; i++) blks[nblk++] = (size_t)atoi(argv[i]) * 1024;
    } else {
        size_t def[] = {64, 256, 1024, 4096, 8192, 65536};
        for (unsigned i = 0; i < sizeof(def) / sizeof(def[0]); i++) blks[nblk++] = def[i] * 1024;
    }

    FILE *f = fopen(argv[1], "rb");
    if (!f) { perror(argv[1]); return 1; }
    fseek(f, 0, SEEK_END);
    long n = ftell(f);
    fseek(f, 0, SEEK_SET);
    uint8_t *src = malloc((size_t)n);
    if (!src || fread(src, 1, (size_t)n, f) != (size_t)n) { fprintf(stderr, "read\n"); return 1; }
    fclose(f);

    int fails = 0;
    for (int b = 0; b < nblk; b++) {
        size_t blk = blks[b];
        size_t bound = cpenc_max_compressed_size(blk);
        uint8_t *comp = malloc(bound), *dec = malloc(blk);
        int bad = 0;
        for (size_t off = 0; off < (size_t)n; off += blk) {
            size_t len = (size_t)n - off < blk ? (size_t)n - off : blk;
            size_t cl = bound, dl = blk;
            if (cpenc_compress(src + off, len, comp, &cl, order) != 0) {
                printf("blk=%zuK off=%zu compress FAILED\n", blk / 1024, off);
                bad++; continue;
            }
            uint8_t stride = ((CPEncHeader *)comp)->stride;
            if (cpenc_decompress(comp, cl, dec, &dl) != 0) {
                printf("blk=%zuK off=%zu len=%zu stride=%s decompress FAILED\n",
                       blk / 1024, off, len, sname(stride));
                bad++;
            } else if (dl != len) {
                printf("blk=%zuK off=%zu len=%zu stride=%s length mismatch: got %zu\n",
                       blk / 1024, off, len, sname(stride), dl);
                bad++;
            } else if (memcmp(dec, src + off, len) != 0) {
                size_t d = 0;
                while (d < len && dec[d] == src[off + d]) d++;
                printf("blk=%zuK off=%zu len=%zu stride=%s BYTES DIFFER first at +%zu "
                       "(got %02x want %02x)\n",
                       blk / 1024, off, len, sname(stride), d, dec[d], src[off + d]);
                bad++;
            }
            if (bad >= 5) { printf("  ...stopping after 5 failures at this block size\n"); break; }
        }
        if (!bad) printf("blk=%zuK  OK (%zu blocks)\n", blk / 1024, ((size_t)n + blk - 1) / blk);
        fails += bad;
        free(comp); free(dec);
    }
    free(src);
    printf("%s: %d failures\n", argv[1], fails);
    return fails ? 1 : 0;
}
