/* fuzz_seed.c — write seeds for fuzz_cpenc.
 *
 * Without seeds the fuzzer spends its budget rediscovering the header from
 * random bytes and never reaches the deep paths.  These samples are chosen
 * to hit each stride the encoder can pick: RAW (incompressible), HEX (hex
 * text), SEG (structured), BWT (prose), REF/REF2 (repetition), plus zeros.
 *
 * Input layout matches the target: [u32 capacity][block bytes].
 *
 *   usage: fuzz_seed <corpus_dir>
 */
#include "../cpenc.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MB (1024 * 1024)

static uint64_t s = 0x123456789ABCDEFull;
static uint32_t rnd(void) { s ^= s << 13; s ^= s >> 7; s ^= s << 17; return (uint32_t)(s >> 32); }

static int seed_no;

static void emit(const uint8_t *src, size_t n, int order, size_t blk, const char *dir) {
    size_t cap = cpenc_max_compressed_size(blk);
    uint8_t *c = malloc(cap);
    if (!c) return;
    for (size_t off = 0; off + 64 <= n; off += blk) {
        size_t len = n - off < blk ? n - off : blk;
        size_t cl = cap;
        if (cpenc_compress(src + off, len, c, &cl, order) != 0) continue;
        char path[512];
        snprintf(path, sizeof(path), "%s/seed%03d_%zu_%d", dir, seed_no++, blk, order);
        FILE *f = fopen(path, "wb");
        if (!f) continue;
        uint32_t room = (uint32_t)len;       /* the capacity a correct caller uses */
        fwrite(&room, 4, 1, f);
        fwrite(c, 1, cl, f);
        fclose(f);
    }
    free(c);
}

int main(int argc, char **argv) {
    const char *dir = argc > 1 ? argv[1] : "/tmp/fz_cpenc";
    size_t N = 256 * 1024;
    uint8_t *b = malloc(N);
    if (!b) return 1;

    /* incompressible */
    for (size_t i = 0; i < N; i++) b[i] = (uint8_t)rnd();
    for (int o = 0; o <= 4; o += 2) { emit(b, N, o, 65536, dir); emit(b, N, o, 4096, dir); }

    /* hex text: "deadbeef" lines */
    for (size_t i = 0; i < N; i += 40) snprintf((char *)b + i, 40, "%08x%08x%08x%08x\n", rnd(), rnd(), rnd(), rnd());
    for (int o = 0; o <= 4; o += 2) emit(b, N, o, 65536, dir);

    /* prose (BWT territory) */
    const char *w = "the quick brown fox jumps over the lazy dog while the ";
    for (size_t i = 0; i < N; i++) b[i] = (uint8_t)w[i % strlen(w)];
    for (int o = 0; o <= 4; o += 2) { emit(b, N, o, 65536, dir); emit(b, N, o, 8192, dir); }

    /* structured numeric records (SEG / REF2 territory) */
    for (size_t i = 0, r = 0; i + 16 <= N; i += 16, r++) {
        int64_t t = 1792300000000LL + (int64_t)r * 1000;
        int32_t id = (int32_t)(r % 50000), q = (int16_t)(r % 300), fl = 7;
        memcpy(b + i, &t, 8); memcpy(b + i + 8, &id, 4);
        memcpy(b + i + 12, &q, 2); memcpy(b + i + 14, &fl, 2);
    }
    for (int o = 0; o <= 4; o += 2) emit(b, N, o, 65536, dir);

    /* highly repetitive (REF / REF2 territory) */
    for (size_t i = 0; i < N; i++) b[i] = (uint8_t)("2026-10-21 [INFO] request id=%08x status=200\n"[i % 45]);
    for (int o = 0; o <= 4; o += 2) emit(b, N, o, 65536, dir);

    /* zeros and long runs (tag 0 / RLE) */
    memset(b, 0, N);
    for (int o = 0; o <= 4; o += 2) emit(b, N, o, 65536, dir);
    for (size_t i = 0; i < N; i++) b[i] = (uint8_t)((i / 500) & 0xFF);
    for (int o = 0; o <= 4; o += 2) emit(b, N, o, 65536, dir);

    printf("wrote %d seeds to %s\n", seed_no, dir);
    free(b);
    return 0;
}
