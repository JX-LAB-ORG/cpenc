/* cpenc_capacity_test — the decoder must never write past the caller's *dl.
 *
 * cpenc_decompress(src, sl, dst, &dl) takes the destination capacity in *dl.
 * A caller sizes that buffer from its own metadata — cfs uses the inode
 * table's orig_size — which a hostile container can make smaller than the
 * block header's original_size.  Every branch must therefore refuse, not
 * write.  The BWT branch did write: memcpy(d, tmp + 1, original_size) with no
 * capacity check, which is a heap overflow reachable from a tampered .cfs.
 *
 * Two halves:
 *   a. honest blocks, dishonest capacity — encode real data with every model
 *      the encoder picks, then decode it into buffers far smaller than the
 *      block declares;
 *   b. dishonest blocks — hand-built streams whose declared sizes do not match
 *      what they decode to, which is what an attacker actually supplies.
 *
 * The buffers sit flush against a PROT_NONE page, the project's stand-in for
 * AddressSanitizer (which cannot run under PRoot).  A write past the capacity
 * faults, the handler longjmps back, and the case is reported as a failure
 * instead of taking the whole test down.
 */
#include "../cpenc.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <setjmp.h>
#include <signal.h>
#include <sys/mman.h>

#define PG 4096

static int failed, checks;
#define CHK(c,m) do{checks++;if(!(c)){failed++;fprintf(stderr,"  FAIL: %s\n",m);}else{fprintf(stderr,"  OK %s\n",m);}}while(0)

static sigjmp_buf fault_jb;
static volatile sig_atomic_t fault_hit;

static void on_fault(int sig) {
    (void)sig;
    fault_hit = 1;
    siglongjmp(fault_jb, 1);
}

static uint32_t rng = 0x12345678;
static uint32_t rnd(void) {
    rng ^= rng << 13; rng ^= rng >> 17; rng ^= rng << 5;
    return rng;
}

/* Buffer whose last byte is flush against an unmapped page. */
static uint8_t *guard_buf(size_t len, size_t *mapped_len, void **map_base) {
    size_t pages = (len + PG - 1) / PG;
    if (pages == 0) pages = 1;
    size_t total = (pages + 2) * PG;
    uint8_t *base = mmap(NULL, total, PROT_READ | PROT_WRITE,
                         MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    if (base == MAP_FAILED) return NULL;
    mprotect(base, PG, PROT_NONE);
    mprotect(base + (pages + 1) * PG, PG, PROT_NONE);
    *mapped_len = total;
    *map_base = base;
    return base + PG + pages * PG - len;
}

static uint8_t *enc(const uint8_t *src, size_t n, int order, size_t *cl) {
    size_t cap = cpenc_max_compressed_size(n), c = cap;
    uint8_t *cb = malloc(cap);
    if (!cb) return NULL;
    if (cpenc_compress(src, n, cb, &c, order) != 0) { free(cb); return NULL; }
    *cl = c;
    return cb;
}

/* Decode into exactly `cap` usable bytes; returns the decoder's rc, or -999
 * if the decoder wrote past the buffer (caught by the guard page). */
static int decode_with_capacity(const uint8_t *cb, size_t cl, size_t cap) {
    size_t ml; void *mb;
    uint8_t *dst = guard_buf(cap ? cap : 1, &ml, &mb);
    if (!dst) return -998;
    size_t dl = cap;
    fault_hit = 0;
    int rc;
    if (sigsetjmp(fault_jb, 1) == 0) {
        rc = cpenc_decompress(cb, cl, dst, &dl);
    } else {
        rc = -999;                      /* the fault handler jumped back here */
    }
    munmap(mb, ml);
    return fault_hit ? -999 : rc;
}

/* Source-like text with a large vocabulary: locally predictable syntax, no
 * long exact repeats — the shape the encoder chooses BWT for.  A periodic
 * sentence would be taken by the reference layer instead. */
static void gen_source_text(uint8_t *buf, size_t n) {
    static const char *syl[] = {
        "ka","to","ri","men","sa","lo","vin","dra","pe","qu",
        "an","el","or","ist","ure","con","ter","pre","dis","tion"
    };
    size_t o = 0;
    while (o < n) {
        char name[64], line[160];
        size_t k = 0;
        for (int w = 0; w < 3 && k < sizeof name - 8; w++) {
            const char *s = syl[rnd() % 20];
            while (*s && k < sizeof name - 1) name[k++] = *s++;
        }
        name[k] = 0;
        for (char *p = name; *p; p++) if (rnd() & 1) *p = (char)(*p - 'a' + 'A');
        int len = snprintf(line, sizeof line, "#define %s %lu\n",
                           name, (unsigned long)(rnd() % 99999));
        size_t take = (size_t)len;
        if (take > n - o) take = n - o;
        memcpy(buf + o, line, take);
        o += take;
    }
}

/* ── honest blocks, dishonest capacity ────────────────────────────────── */
static int covered_bwt;

static void try_dataset(const char *name, const uint8_t *data, size_t n) {
    size_t cl;
    uint8_t *cb = enc(data, n, 2, &cl);
    if (!cb) { failed++; checks++; fprintf(stderr, "  FAIL: %s: compress failed\n", name); return; }
    unsigned stride = cb[3];
    printf("  [%s] n=%zu stride=0x%02X cl=%zu\n", name, n, stride, cl);
    if (stride == CPENC_STRIDE_BWT) covered_bwt = 1;

    size_t caps[] = { 0, 1, 16, n / 2, n > 1 ? n - 1 : 0 };
    for (size_t i = 0; i < sizeof caps / sizeof caps[0]; i++) {
        if (caps[i] >= n) continue;
        char what[128];
        snprintf(what, sizeof what, "%s: capacity %zu < %zu must be refused", name, caps[i], n);
        int rc = decode_with_capacity(cb, cl, caps[i]);
        checks++;
        if (rc == -999) { failed++; fprintf(stderr, "  FAIL: %s — wrote past the capacity\n", what); }
        else if (rc == 0) { failed++; fprintf(stderr, "  FAIL: %s — returned success\n", what); }
        else fprintf(stderr, "  OK %s\n", what);
    }

    /* exact capacity: must succeed and reproduce the input */
    size_t ml; void *mb;
    uint8_t *dst = guard_buf(n ? n : 1, &ml, &mb);
    if (!dst) { failed++; checks++; free(cb); return; }
    size_t dl = n;
    int rc = cpenc_decompress(cb, cl, dst, &dl);
    checks++;
    if (rc != 0 || dl != n || memcmp(dst, data, n)) {
        failed++;
        fprintf(stderr, "  FAIL: %s: exact-capacity decode (rc=%d dl=%zu)\n", name, rc, dl);
    }
    munmap(mb, ml);
    free(cb);
}

/* ── dishonest blocks ─────────────────────────────────────────────────── */
/* A hex wrapper whose declared text length is tiny while its packed buffer is
 * large: the restore loop would write 2*bin_len bytes into a 4-byte buffer. */
static void dishonest_hex(void) {
    enum { BIN = 4096, RESTORE = 4 };
    size_t paylen = 13 + BIN;
    uint8_t *pay = calloc(1, paylen);
    if (!pay) { failed++; checks++; return; }
    pay[0] = RESTORE;                       /* orig_len, little-endian */
    pay[5] = BIN & 0xFF; pay[6] = (BIN >> 8) & 0xFF;   /* bin_len, n_sep = 0 */
    size_t icl;
    uint8_t *inner = enc(pay, paylen, 0, &icl);
    if (!inner) { free(pay); failed++; checks++; return; }
    unsigned rs = inner[3];
    if (rs == CPENC_STRIDE_HEX || rs == CPENC_STRIDE_SEG) rs = CPENC_STRIDE_RAW;

    size_t ol = CPENC_HEADER_SIZE + 1 + icl;
    uint8_t *outer = calloc(1, ol);
    uint8_t *op = outer;
    *op++ = (uint8_t)(CPENC_MAGIC & 0xFF); *op++ = (uint8_t)(CPENC_MAGIC >> 8);
    *op++ = CPENC_VERSION;
    *op++ = CPENC_STRIDE_HEX;
    uint32_t o_len = (uint32_t)paylen, o_comp = (uint32_t)(1 + icl);
    *op++ = o_len & 0xFF; *op++ = (o_len >> 8) & 0xFF; *op++ = (o_len >> 16) & 0xFF; *op++ = (o_len >> 24) & 0xFF;
    *op++ = o_comp & 0xFF; *op++ = (o_comp >> 8) & 0xFF; *op++ = (o_comp >> 16) & 0xFF; *op++ = (o_comp >> 24) & 0xFF;
    *op++ = (uint8_t)rs;
    memcpy(op, inner, icl);

    int rc = decode_with_capacity(outer, ol, RESTORE);
    checks++;
    if (rc == -999) { failed++; fprintf(stderr, "  FAIL: dishonest hex wrote past the capacity\n"); }
    else if (rc == 0) { failed++; fprintf(stderr, "  FAIL: dishonest hex returned success\n"); }
    else fprintf(stderr, "  OK dishonest hex rejected\n");

    free(outer); free(inner); free(pay);
}

/* A BWT block whose header claims N bytes while the inner stream decodes to
 * N+1: the unchecked memcpy would write N bytes from a 4-byte destination. */
static void dishonest_bwt(void) {
    enum { N = 32768 };
    uint8_t *raw = calloc(1, N + 1);       /* the bytes the inner stream yields */
    if (!raw) { failed++; checks++; return; }
    size_t icl;
    uint8_t *inner = enc(raw, N + 1, 0, &icl);
    if (!inner) { free(raw); failed++; checks++; return; }

    size_t ol = CPENC_HEADER_SIZE + 4 + icl;
    uint8_t *outer = calloc(1, ol);
    uint8_t *op = outer;
    *op++ = (uint8_t)(CPENC_MAGIC & 0xFF); *op++ = (uint8_t)(CPENC_MAGIC >> 8);
    *op++ = CPENC_VERSION;
    *op++ = CPENC_STRIDE_BWT;
    uint32_t o_orig = N, o_comp = (uint32_t)(4 + icl);
    *op++ = o_orig & 0xFF; *op++ = (o_orig >> 8) & 0xFF; *op++ = (o_orig >> 16) & 0xFF; *op++ = (o_orig >> 24) & 0xFF;
    *op++ = o_comp & 0xFF; *op++ = (o_comp >> 8) & 0xFF; *op++ = (o_comp >> 16) & 0xFF; *op++ = (o_comp >> 24) & 0xFF;
    *op++ = 0; *op++ = 0; *op++ = 0; *op++ = 0;          /* primary index */
    memcpy(op, inner, icl);

    int rc = decode_with_capacity(outer, ol, 4);
    checks++;
    if (rc == -999) { failed++; fprintf(stderr, "  FAIL: dishonest BWT wrote past the capacity\n"); }
    else if (rc == 0) { failed++; fprintf(stderr, "  FAIL: dishonest BWT returned success\n"); }
    else fprintf(stderr, "  OK dishonest BWT rejected\n");

    free(outer); free(inner); free(raw);
}

int main(void) {
    struct sigaction sa;
    memset(&sa, 0, sizeof sa);
    sa.sa_handler = on_fault;
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = SA_NODEFER;
    sigaction(SIGSEGV, &sa, NULL);
    sigaction(SIGBUS, &sa, NULL);

    printf("[cap1] honest blocks, capacity smaller than declared\n");
    {
        size_t n = 65536;
        uint8_t *defs = malloc(n);
        gen_source_text(defs, n);
        try_dataset("source-like text", defs, n);
        free(defs);
    }
    {
        size_t n = 131072;
        uint8_t *hex = malloc(n);
        for (size_t i = 0; i < n; i++) {
            static const char h[] = "0123456789ABCDEF";
            hex[i] = (i % 33 == 32) ? '\n' : (uint8_t)h[rnd() & 15];
        }
        try_dataset("uppercase hex", hex, n);
        free(hex);
    }
    {
        size_t n = 65536;
        uint8_t *rn = malloc(n);
        for (size_t i = 0; i < n; i++) rn[i] = (uint8_t)(rnd() & 0xFF);
        try_dataset("random", rn, n);
        free(rn);
    }
    {
        size_t n = 262144;
        uint8_t *rec = malloc(n);
        for (size_t i = 0; i + 32 <= n; i += 32) {
            uint32_t k = rnd() % 64;
            memcpy(rec + i, &k, 4); memcpy(rec + i + 4, &k, 4);
            memset(rec + i + 8, 0, 24);
        }
        try_dataset("records", rec, n);
        free(rec);
    }
    {
        /* two megabytes of mixed content: the shape that makes the encoder
         * segment, which is the branch with per-segment buffers */
        size_t n = 2u << 20;
        uint8_t *mix = malloc(n);
        for (size_t i = 0; i < n; i++)
            mix[i] = (i < n / 2) ? (uint8_t)(rnd() & 0xFF)
                                 : (uint8_t)("structured text with repeats, repeats, repeats\n"[i % 47]);
        try_dataset("mixed 2MB", mix, n);
        free(mix);
    }
    CHK(covered_bwt, "the BWT branch was exercised (the regression this file exists for)");

    printf("[cap2] dishonest blocks\n");
    dishonest_hex();
    dishonest_bwt();
    {
        /* The fuzzer's find: a valid BWT block whose 4-byte primary index is
         * replaced by a value far outside the block.  bwt_dec walked the LF
         * chain from it and read past the block. */
        size_t n = 32768;
        uint8_t *src = malloc(n);
        gen_source_text(src, n);
        size_t cap_len;
        uint8_t *cb = enc(src, n, 2, &cap_len);
        checks++;
        if (!cb || cb[3] != CPENC_STRIDE_BWT) {
            failed++;
            fprintf(stderr, "  FAIL: dishonest BWT index: encoder did not pick BWT\n");
        } else {
            cb[CPENC_HEADER_SIZE + 0] = 0x09;
            cb[CPENC_HEADER_SIZE + 1] = 0x15;
            cb[CPENC_HEADER_SIZE + 2] = 0xff;
            cb[CPENC_HEADER_SIZE + 3] = 0xff;
            int rc = decode_with_capacity(cb, cap_len, n);
            if (rc == -999) { failed++; fprintf(stderr, "  FAIL: bad BWT index wrote past the capacity\n"); }
            else if (rc == 0) { failed++; fprintf(stderr, "  FAIL: bad BWT index returned success\n"); }
            else fprintf(stderr, "  OK bad BWT index rejected\n");
        }
        free(cb); free(src);
    }

    printf("\n=== %d failures ===\n", failed);
    return failed ? 1 : 0;
}
