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

/* A BWT block paired with a capacity far smaller than the block declares.
 *
 * The BWT branch rebuilds h->original_size bytes from the inner stream; a
 * caller that sizes its output buffer from some other (equally untrusted)
 * field — the inode table is exactly such a place — gets an overflow unless
 * the decoder compares the two.  Random mutation cannot make this input:
 * the block has to be a genuinely self-consistent BWT stream, which only the
 * encoder can produce, and then have its capacity prefix disagree with it.
 *
 * The source text is synthetic source code on purpose.  Prose and repeated
 * logs go to REF/REF2, and structured tables go to the delta path; code-like
 * text with a wide identifier vocabulary is what actually makes the encoder
 * pick stride 0xFC.  Verified below before anything is written.
 */
#define BWT_SEED_N  16384
#define BWT_SEED_CAP 16u

static void emit_bwt_small_cap(const char *dir) {
    static const char *kw[] = {
        "int", "return", "static", "if", "else", "for", "while", "struct",
        "const", "char", "size_t", "void", "uint32_t", "buf", "len", "data",
        "ptr", "idx", "count", "value"
    };
    uint8_t *src = malloc(BWT_SEED_N);
    size_t cap = cpenc_max_compressed_size(BWT_SEED_N);
    uint8_t *c = malloc(cap);
    if (!src || !c) { free(src); free(c); return; }

    /* deterministic, independent of the rnd() stream used above */
    uint64_t rs = 0xD1B54A32D192ED03ull;
    size_t o = 0;
    while (o < BWT_SEED_N) {
        size_t room = BWT_SEED_N - o;
        const char *a = kw[rs % 20], *b2 = kw[(rs >> 8) % 20], *e = kw[(rs >> 16) % 20];
        unsigned x = (unsigned)(rs % 1000), y = (unsigned)((rs >> 20) % 100000);
        int w;
        rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17;
        switch (rs % 5) {
        case 0: w = snprintf((char *)src + o, room, "    %s %s_%u = %u;\n", a, b2, x, y); break;
        case 1: w = snprintf((char *)src + o, room, "if (%s_%u > %u) { %s(%s, %u); }\n", a, x % 100, y, b2, e, x); break;
        case 2: w = snprintf((char *)src + o, room, "/* %u %u %u %u */\n", x % 10, y % 10, (x + y) % 10, (x * y) % 10); break;
        case 3: w = snprintf((char *)src + o, room, "%s_%u += %s[%u] * 0x%x;\n", a, x % 50, b2, y % 64, x % 256); break;
        default:w = snprintf((char *)src + o, room, "while (%s < %u) { %s_%u++; }\n", a, x, b2, y % 77); break;
        }
        if (w < 0) break;
        o += (size_t)w < room ? (size_t)w : room;   /* snprintf truncates at room */
    }

    size_t cl = cap;
    CPEncHeader h;
    if (o != BWT_SEED_N || cpenc_compress(src, BWT_SEED_N, c, &cl, 2) != 0) {
        fprintf(stderr, "bwt-small-cap: encoder call failed\n");
        free(src); free(c);
        return;
    }
    memcpy(&h, c, sizeof h);
    if (h.stride != CPENC_STRIDE_BWT && h.stride != CPENC_STRIDE_BWT2) {
        /* Force-BWT failed: say so rather than writing a seed that misses
         * the branch.  (Not expected with this text; order 4 is the next
         * knob if a future encoder gate rejects order 2.) */
        fprintf(stderr, "bwt-small-cap: encoder picked stride 0x%02X, not a BWT stride; seed not written\n",
                h.stride);
        free(src); free(c);
        return;
    }

    char path[512];
    snprintf(path, sizeof path, "%s/bwt-small-cap.seed", dir);
    FILE *f = fopen(path, "wb");
    if (!f) { free(src); free(c); return; }
    uint32_t room = BWT_SEED_CAP;               /* deliberately too small */
    fwrite(&room, 4, 1, f);
    fwrite(c, 1, cl, f);
    fclose(f);
    printf("wrote %s (%zu bytes, stride 0x%02X, declares %u)\n",
           path, cl + 4, h.stride, h.original_size);
    free(src); free(c);
}

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

    emit_bwt_small_cap(dir);

    printf("wrote %d seeds to %s\n", seed_no, dir);
    free(b);
    return 0;
}
