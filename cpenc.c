#include "cpenc.h"
#include <pthread.h>
#include "bwt.h"
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <stdio.h>

#if defined(__aarch64__)
#include "cpenc_neon.h"
#endif
#if defined(__x86_64__) && defined(__AVX2__)
#include "cpenc_avx2.h"
#endif


/* ── Header writes ───────────────────────────────────────────────
 * The destination is the caller's buffer, and an internal nested block is
 * not guaranteed to start on a 4-byte boundary: a segmented block written
 * inside a hex wrapper begins right after the wrapper's one-byte stride
 * marker.  Writing CPEncHeader fields through a struct pointer there is an
 * unaligned store — undefined behaviour (UBSan flags it, strict-alignment
 * targets fault).  These stay byte-wise, mirroring what the decoder already
 * does with memcpy on the read side. */
static inline void hdr_head(uint8_t *p, uint8_t stride, uint32_t orig) {
    p[0] = (uint8_t)(CPENC_MAGIC & 0xFF);
    p[1] = (uint8_t)(CPENC_MAGIC >> 8);
    p[2] = (uint8_t)CPENC_VERSION;
    p[3] = stride;
    p[4] = (uint8_t)orig; p[5] = (uint8_t)(orig >> 8);
    p[6] = (uint8_t)(orig >> 16); p[7] = (uint8_t)(orig >> 24);
}
static inline void hdr_comp(uint8_t *p, uint32_t comp) {
    p[8]  = (uint8_t)comp; p[9]  = (uint8_t)(comp >> 8);
    p[10] = (uint8_t)(comp >> 16); p[11] = (uint8_t)(comp >> 24);
}
static inline void hdr_put(uint8_t *p, uint8_t stride, uint32_t orig, uint32_t comp) {
    hdr_head(p, stride, orig);
    hdr_comp(p, comp);
}

#define HUFF_MAX 12
#define HUFF_WSZ (1<<HUFF_MAX)

/* Hard ceiling on a block's decompressed size.  The format's field is a
 * u32, but nothing in cpenc or cfs ever produces a block anywhere near this:
 * cfs splits files at 64MB and the device layer uses 8MB.  Without a ceiling
 * a 28-byte block claiming 1GB makes the decoder allocate and grind for
 * seconds — fuzzing found that, and it is a denial-of-service knob that a
 * container controls. */
#define CPENC_MAX_BLOCK (256u << 20)

/* How many bytes are left between a cursor and the end of its buffer.
 * `(size_t)(end - p)` is the usual way to write this and it is wrong once p
 * is already past end: the subtraction wraps, the comparison passes, and the
 * "not enough bytes left" guard lets a read through.  Fuzzing found exactly
 * that — h_dec was handed an end 256 bytes BEHIND its cursor and read into a
 * guard page.  Every length check in the decoder uses this now. */
#define REMAIN(p, end) ((p) > (end) ? (size_t)0 : (size_t)((end) - (p)))
#define HUFF_TABSZ (HUFF_WSZ * 2)

/* ── entropy from freq table ─────────────────────────────────
 * fq_ent is called per group (tens of thousands per block); 256 libm
 * log2 calls each would dominate — precompute log2 for the integer
 * count range (segmented groups are <= 64KB so the table covers the
 * hot path; larger counts fall back to log2). */
static double lg2t[65537];
static int lg2_init_done = 0;
static void lg2_init(void) {
    if (lg2_init_done) return;
    for (int i = 0; i <= 65536; i++) lg2t[i] = log2((double)i);
    lg2_init_done = 1;
}
static double fq_ent(const uint32_t fq[256], size_t tot) {
    if(!tot) return 0;
    double h=0;
    double lgt = tot <= 65536 ? lg2t[tot] : log2((double)tot);
    for(int i=0;i<256;i++){
        uint32_t c = fq[i];
        if(c){
            double lc = c <= 65536 ? lg2t[c] : log2((double)c);
            h -= (double)c * (lc - lgt) / (double)tot;
        }
    }
    return h;
}

static double ent(const uint8_t *p, size_t n) {
    uint32_t c[256]={0};
    size_t k = n < 8192 ? n : 8192;
    for(size_t i=0;i<k;i++) c[p[i]]++;
    return fq_ent(c, k);
}

static void tr(const uint8_t *s, uint8_t *d, size_t n, int str) {
    size_t e=n/(size_t)str, end=e*(size_t)str;
    for(int g=0;g<str;g++)
        for(size_t i=(size_t)g;i<end;i+=(size_t)str) *d++=s[i];
    if(end<n) memcpy(d, s+end, n-end);
}
static void ur(const uint8_t *s, uint8_t *d, size_t n, int str) {
    size_t e=n/(size_t)str, end=e*(size_t)str;
    for(int g=0;g<str;g++){
        size_t o=(size_t)g*e;
        for(size_t i=(size_t)g;i<end;i+=(size_t)str) d[i]=s[o++];
    }
    if(end<n) memcpy(d+end, s+end, n-end);
}

static void de(uint8_t *p, size_t n, int o) {
    if(n<2) return;
    while(o--) for(size_t i=n-1;i>0;i--) p[i]-=p[i-1];
}
static void dd(uint8_t *p, size_t n, int o) {
    if(n<2) return;
#if defined(__aarch64__) && !defined(CPENC_SCALAR_DD)
    /* Prefix sum via path doubling inside each 16-byte lane (Kogge-Stone):
     * four shift+add steps replace the 16-step serial carry chain. */
    while(o--) ne_undelta1(p, n);
#else
    while(o--) for(size_t i=1;i<n;i++) p[i]+=p[i-1];
#endif
}

static int ds(const uint8_t *p, size_t n) {
    size_t k=n<1024?n:1024;
    uint8_t t[1024], d[1024];
    int bs=1; double be=1e99;
    int cd[4]={1,2,4,8};
    for(int ci=0;ci<4;ci++){
        int s=cd[ci]; if(k<(size_t)s*4) continue;
        size_t e=k/(size_t)s;
        tr(p,t,k,s);
        for(int g=0;g<s;g++){uint8_t *q=t+(size_t)g*e; d[g*e]=q[0];
            for(size_t j=1;j<e;j++) d[g*e+j]=q[j]-q[j-1];}
        double ev=ent(d,k); if(ev<be){be=ev;bs=s;}
    }
    return bs;
}

/* ── vector-space features (algebraic classification) ─────────────
 * The block is a vector x over Z_256.  Instead of trial-compressing every
 * candidate, measure the quantities the linear operators act on:
 *
 *   z    = #{i : x[i] == x[i-1]} / (n-1)      sparsity of the difference
 *                                              vector Dx  -> delta family
 *   a_k  = #{i : x[i] == x[i-k]} / (n-k)      lag-k equality (autocorrelation
 *                                              of the value vector) -> stride
 *   w    = repeated 4-byte windows            collision rate of the window
 *                                              vector w_i -> match layer
 *
 * w is exactly the "f(x1)=f(x2)" statistic: it counts how often two
 * positions of the window vector take the same value.  All of this is
 * O(n) counting on a bounded sample, with no histogram and no trial. */
typedef struct { double z, a1, a2, a4, a8, w, lr, lr16, lr32, cover, ap1, ap2, rsim, rdup, dapprox, dshape, ml_short, ml_long; double rho[17]; int kbest; double drho; } VFeat;

/* with_rho: the lag-k correlation is O(16n) 64-bit multiplies and must not
 * run on the hot path; the match-layer gate needs only w, which costs one
 * pass over a 64KB window sample. */
static void vfeat(const uint8_t *p, size_t n, VFeat *f, int with_rho) {
    size_t k = n < (64u << 10) ? n : (64u << 10);
    memset(f, 0, sizeof *f);
    if (k < 16) return;
    size_t z = 0, a2 = 0, a4 = 0, a8 = 0;
    for (size_t i = 1; i < k; i++) if (p[i] == p[i-1]) z++;
    for (size_t i = 2; i < k; i++) if (p[i] == p[i-2]) a2++;
    for (size_t i = 4; i < k; i++) if (p[i] == p[i-4]) a4++;
    for (size_t i = 8; i < k; i++) if (p[i] == p[i-8]) a8++;
    f->z = (double)z / (double)(k - 1);
    f->a1 = f->z;
    f->a2 = (double)a2 / (double)(k - 2);
    f->a4 = (double)a4 / (double)(k - 4);
    f->a8 = (double)a8 / (double)(k - 8);

    /* Correlation of the value vector with its lag-k shift.  The peak of
     * rho[k]-rho[k-1] is the structural period (stride); the decay rate is
     * how predictable the block is.  Measured at 512K: f32 peaks at k=4
     * (d=0.17), ELF at k=8 (d=0.10), text is flat. */
    if (with_rho) {
        size_t m = n < (16u << 10) ? n : (16u << 10);
        for (int lag = 1; lag <= 8; lag++) {
            if ((size_t)lag + 16 > m) break;
            uint64_t sxy = 0, sxx = 0, syy = 0;
            for (size_t i = (size_t)lag; i < m; i++) {
                uint64_t a = p[i], b = p[i - lag];
                sxy += a * b; sxx += a * a; syy += b * b;
            }
            double den = sqrt((double)sxx) * sqrt((double)syy);
            f->rho[lag] = den > 0 ? (double)sxy / den : 0;
        }
        for (int lag = 2; lag <= 8; lag++) {
            double d = f->rho[lag] - f->rho[lag - 1];
            if (d > 0.02 && d > f->drho) { f->drho = d; f->kbest = lag; }
        }
    }

    enum { WB = 16, WS = 1u << WB };
    uint64_t *tab = calloc(WS, sizeof(uint64_t));
    size_t rep = 0, tot = 0;
    if (tab) {
        for (size_t i = 0; i + 4 <= k; i += 4) {
            uint32_t v; memcpy(&v, p + i, 4);
            uint64_t key = ((uint64_t)v << 1) | 1u;
            uint32_t h = (v * 2654435761u) >> (32 - WB);
            for (int pr = 0; pr < 8; pr++) {
                uint64_t *slot = &tab[(h + (uint32_t)pr) & (WS - 1)];
                if (*slot == 0) { *slot = key; break; }
                if (*slot == key) { rep++; break; }
            }
            tot++;
        }
        free(tab);
    }
    f->w = tot ? (double)rep / (double)tot : 0;

    /* Region-pair statistics: split the sample into 1KB regions, describe
     * each by (mean value, net slope) and a hash.  rdup = fraction of pairs
     * that are byte-identical (REF's territory); rsim = fraction whose mean
     * and slope are both close (BWT's territory: two different regions with
     * approximately the same f and k). */
    {
        enum { RG = 1024, RM = 64 };
        size_t nr = k / RG; if (nr > RM) nr = RM;
        if (nr >= 2) {
            double mu[RM]; double sl[RM]; uint64_t hh[RM];
            for (size_t r = 0; r < nr; r++) {
                const uint8_t *q = p + r * RG;
                uint64_t sum = 0, h = 0xcbf29ce484222325ULL;
                for (size_t i = 0; i < RG; i++) {
                    sum += q[i];
                    h = (h ^ q[i]) * 0x100000001b3ULL;
                }
                int net = (int)q[RG-1] - (int)q[0];
                if (net > 127) net -= 256; if (net < -128) net += 256;
                mu[r] = (double)sum / RG;
                sl[r] = (double)net / RG;
                hh[r] = h;
            }
            size_t dup = 0, sim = 0, pairs = 0;
            for (size_t a = 0; a < nr; a++)
                for (size_t b = a + 1; b < nr; b++) {
                    pairs++;
                    if (hh[a] == hh[b]) dup++;
                    double dm = mu[a] - mu[b]; if (dm < 0) dm = -dm;
                    double ds = sl[a] - sl[b]; if (ds < 0) ds = -ds;
                    if (dm <= 8 && ds <= 0.05) sim++;
                }
            if (pairs) { f->rdup = (double)dup / pairs; f->rsim = (double)sim / pairs; }
        }
    }

    if (with_rho) {
    /* ap2: curvature — how often the slope itself barely changes
     * (|k_i - k_{i-1}| <= 2).  BWT needs BOTH this and ap1: a nearly flat
     * level set that also varies slowly is what turns into long MTF runs. */
    {
        size_t c2 = 0;
        for (size_t i = 2; i < k; i++) {
            int k1 = (int)p[i] - (int)p[i-1];
            int k0 = (int)p[i-1] - (int)p[i-2];
            if (k1 > 127) k1 -= 256; if (k1 < -128) k1 += 256;
            if (k0 > 127) k0 -= 256; if (k0 < -128) k0 += 256;
            int d = k1 - k0; if (d < 0) d = -d;
            if (d <= 2) c2++;
        }
        f->ap2 = (double)c2 / (double)(k - 2);
    }

    /* ap1: how often adjacent values are within +-2 (mod 256) — the
     * "k-neighbourhood approximation" that BWT turns into local runs. */
    {
        size_t ap = 0;
        for (size_t i = 1; i < k; i++) {
            int d = (int)p[i] - (int)p[i-1];
            if (d < 0) d = -d;
            if (d > 128) d = 256 - d;
            if (d <= 2) ap++;
        }
        f->ap1 = (double)ap / (double)(k - 1);
    }
    }

    if (with_rho) {
    /* Match-length histogram from a cheap sampled matcher: the direct
     * "short exact vs long exact" quantity.  Greedy over a 64KB sample,
     * 4-byte hash, chain depth 8, match capped at 256 bytes.  Short means
     * 4 .. 15, long means 16+. */
    {
        enum { MH = 1 << 16, MD = 8, MCAP = 256 };
        int32_t *head = malloc(MH * 4), *chain = malloc(MH * 4);
        if (head && chain) {
            memset(head, 0xFF, MH * 4);
            size_t m = k < (64u << 10) ? k : (64u << 10);
            size_t cshort = 0, clong = 0, tot = 0;
            for (size_t i = 0; i + 4 <= m; i++) {
                uint32_t v; memcpy(&v, p + i, 4);
                uint32_t h = ((v * 2654435761u) >> 16) & (MH - 1);
                int32_t c = head[h]; int best = 0, tries = 0;
                while (c >= 0 && tries < MD) {
                    size_t cd = i - (size_t)c;
                    if (cd && !memcmp(p + c, p + i, 4)) {
                        size_t L = 4;
                        while (L < MCAP && i + L < m && cd + L < i && p[c + L] == p[i + L]) L++;
                        if ((int)L > best) best = (int)L;
                    }
                    c = chain[c & (MH - 1)]; tries++;
                }
                chain[i & (MH - 1)] = head[h]; head[h] = (int32_t)i;
                if (best >= 4) { if (best < 16) cshort++; else clong++; }
                tot++;
            }
            f->ml_short = tot ? (double)cshort / (double)tot : 0;
            f->ml_long  = tot ? (double)clong / (double)tot : 0;
        }
        free(head); free(chain);
    }

    /* Similar-region coverage.  Two different x-intervals are "the same
     * region" when f(x) (level), the local slope amplitude |k|, and the
     * amplitude of k's rate of change |k'| agree.  Each 8-byte block is
     * reduced to that signature; repeated signatures at a different position
     * extended over consecutive blocks give the repeat length in x. */
    {
        enum { BS = 8, MAXB = 8192 };
        size_t nb = k / BS; if (nb > MAXB) nb = MAXB;
        uint8_t sig[MAXB], sig2[MAXB];
        for (size_t b = 0; b < nb; b++) {
            const uint8_t *q = p + b * BS;
            unsigned sum = 0, kacc = 0, k2acc = 0;
            for (int j = 0; j < BS; j++) sum += q[j];
            for (int j = 0; j + 1 < BS; j++) {
                int d = (int)q[j+1] - (int)q[j];
                if (d > 127) d -= 256; if (d < -128) d += 256;
                kacc += (unsigned)(d < 0 ? -d : d);
            }
            for (int j = 0; j + 2 < BS; j++) {
                int d1 = (int)q[j+1] - (int)q[j];
                int d2 = (int)q[j+2] - (int)q[j+1];
                if (d1 > 127) d1 -= 256; if (d1 < -128) d1 += 256;
                if (d2 > 127) d2 -= 256; if (d2 < -128) d2 += 256;
                int dd = d2 - d1;
                k2acc += (unsigned)(dd < 0 ? -dd : dd);
            }
            unsigned m1 = kacc / (BS - 1), m2 = k2acc / (BS - 2);
            unsigned fq  = (sum / BS) >> 4;
            unsigned kq  = m1 == 0 ? 0 : m1 < 4 ? 1 : m1 < 16 ? 2 : 3;
            unsigned k2q = m2 == 0 ? 0 : m2 < 4 ? 1 : m2 < 16 ? 2 : 3;
            sig[b] = (uint8_t)(fq * 16 + kq * 4 + k2q);
        }
        int32_t tab[256]; memset(tab, 0xFF, sizeof tab);
        size_t covered = 0, b = 0;
        while (b < nb) {
            int32_t prev = tab[sig[b]];
            tab[sig[b]] = (int32_t)b;
            if (prev >= 0 && (size_t)prev + 8 <= b) {
                size_t d = b - (size_t)prev, len = 1;
                while (b + len < nb && sig[b + len] == sig[b + len - d]) len++;
                if (len >= 8) { covered += len * BS; b += len; continue; }
            }
            b++;
        }
        f->dapprox = nb ? (double)covered / (double)(nb * BS) : 0;

        /* Level-free variant: same derivative profile at a different level is
         * exactly the "base + constant offset" candidate pool. */
        int32_t tab2[256]; memset(tab2, 0xFF, sizeof tab2);
        size_t cov2 = 0; b = 0;
        while (b < nb) {
            int32_t prev = tab2[sig2[b]];
            tab2[sig2[b]] = (int32_t)b;
            if (prev >= 0 && (size_t)prev + 8 <= b) {
                size_t d = b - (size_t)prev, len = 1;
                while (b + len < nb && sig2[b + len] == sig2[b + len - d]) len++;
                if (len >= 8) { cov2 += len * BS; b += len; continue; }
            }
            b++;
        }
        f->dshape = nb ? (double)cov2 / (double)(nb * BS) : 0;
    }

    /* cover: fraction of the sample covered by greedy exact matches of at
     * least 8 bytes — the "same f, different x" quantity REF wins on.
     * O(n) with one hash table; no output is produced. */
    {
        enum { CB = 16, CS = 1u << CB, NIL = 0xFFFFFFFFu };
        uint32_t *pos = malloc(CS * sizeof(uint32_t));
        if (pos) {
            memset(pos, 0xFF, CS * sizeof(uint32_t));
            size_t i = 0, covered = 0;
            while (i + 8 <= k) {
                uint64_t v; memcpy(&v, p + i, 8);
                uint64_t hh = v * 0xff51afd7ed558ccdULL; hh ^= hh >> 33;
                uint32_t h = (uint32_t)hh & (CS - 1);
                uint32_t prev = pos[h];
                pos[h] = (uint32_t)i;
                size_t len = 0;
                if (prev != NIL && prev < i && memcmp(p + prev, p + i, 8) == 0) {
                    while (i + len < k && p[prev + len] == p[i + len]) len++;
                }
                if (len >= 8) { covered += len; i += len; }
                else i++;
            }
            f->cover = (double)covered / (double)k;
            free(pos);
        }
    }
    }

    if (with_rho) {
    /* Long-range repeat rate: 8-byte windows seen before at a distance
     * greater than 4KB.  This is the statistic that separates the blocks
     * REF wins (tables, repeated templates) from prose, where BWT wins. */
    {
        enum { LB = 17, LS = 1u << LB, NIL = 0xFFFFFFFFu };
        uint32_t *pos = malloc(LS * sizeof(uint32_t));
        if (pos) {
            memset(pos, 0xFF, LS * sizeof(uint32_t));
            size_t lrn = 0, lrt = 0;
            for (size_t i = 0; i + 8 <= k; i += 8) {
                uint64_t v; memcpy(&v, p + i, 8);
                v ^= v >> 33; v *= 0xff51afd7ed558ccdULL; v ^= v >> 33;
                uint32_t h = (uint32_t)v & (LS - 1);
                uint32_t prev = pos[h];
                if (prev != NIL && i - prev > 4096) lrn++;
                pos[h] = (uint32_t)i;
                lrt++;
            }
            f->lr = lrt ? (double)lrn / (double)lrt : 0;
            free(pos);
        }
        /* Longer windows: a repeat of 16 or 32 bytes is a long match, which
         * is what REF wins on, as opposed to the short clustering BWT eats. */
        for (int wlen = 16; wlen <= 32; wlen += 16) {
            enum { LS2 = 1u << 17 };
            uint32_t *p2 = malloc(LS2 * sizeof(uint32_t));
            if (!p2) continue;
            memset(p2, 0xFF, LS2 * sizeof(uint32_t));
            size_t rn2 = 0, rt2 = 0;
            for (size_t i = 0; i + (size_t)wlen <= k; i += 8) {
                uint64_t h = 0xcbf29ce484222325ULL;
                for (int b = 0; b < wlen; b++) h = (h ^ p[i + b]) * 0x100000001b3ULL;
                uint32_t idx = (uint32_t)(h ^ (h >> 32)) & (LS2 - 1);
                uint32_t prev = p2[idx];
                if (prev != 0xFFFFFFFFu && i - prev > 4096) rn2++;
                p2[idx] = (uint32_t)i;
                rt2++;
            }
            double v = rt2 ? (double)rn2 / (double)rt2 : 0;
            if (wlen == 16) f->lr16 = v; else f->lr32 = v;
            free(p2);
        }
    }
    }
}

/* 4-interleaved frequency counter — breaks the per-iteration data dependency
 * chain on fq[], letting modern CPUs execute 4 increments in parallel.
 * Measured ~2-3× faster than single-stream fq[p[i]]++ on aarch64. */
static inline __attribute__((always_inline)) void freq256(const uint8_t *p, size_t n, uint32_t fq[256]) {
    uint32_t a[256]={0},b[256]={0},c[256]={0},d[256]={0};
    size_t i=0;
    /* Main loop: 16 bytes per iteration, 4 independent streams */
    for(;i+16<=n;i+=16){
        a[p[i  ]]++; b[p[i+1]]++; c[p[i+2]]++; d[p[i+3]]++;
        a[p[i+4]]++; b[p[i+5]]++; c[p[i+6]]++; d[p[i+7]]++;
        a[p[i+8]]++; b[p[i+9]]++; c[p[i+10]]++; d[p[i+11]]++;
        a[p[i+12]]++; b[p[i+13]]++; c[p[i+14]]++; d[p[i+15]]++;
    }
    for(;i+4<=n;i+=4){a[p[i]]++;b[p[i+1]]++;c[p[i+2]]++;d[p[i+3]]++;}
    /* Tail bytes must land in a[] — writing fq[] directly would be
     * clobbered by the merge below, silently dropping the last 1-3
     * bytes' frequencies (breaks Huffman/LUT for non-4-aligned sizes). */
    for(;i<n;i++) a[p[i]]++;
    /* Merge */
    for(int j=0;j<256;j++) fq[j]=a[j]+b[j]+c[j]+d[j];
}

/* ── Huffman engine ──────────────────────────────────────────── */
typedef struct { uint32_t f; int sym; } HPair;
static void hp_push(HPair *h,int *n,uint32_t f,int sym){
    int i=(*n)++;while(i>0&&h[(i-1)/2].f>f){h[i]=h[(i-1)/2];i=(i-1)/2;}
    h[i].f=f;h[i].sym=sym;
}
static HPair hp_pop(HPair *h,int *n){
    HPair r=h[0];h[0]=h[--(*n)];int i=0;
    while(2*i+1<*n){int c=2*i+1;if(c+1<*n&&h[c+1].f<h[c].f)c++;if(h[i].f<=h[c].f)break;HPair t=h[i];h[i]=h[c];h[c]=t;i=c;}
    return r;
}
static void hlens(const uint32_t fq[256],uint8_t lens[256]){
    memset(lens,0,256);
    uint32_t f[256];
    memcpy(f, fq, sizeof(f));
    /* Code lengths are the tree depths, which are already optimal for the
     * frequencies they were built from — no canonical re-assignment is
     * needed (and one would only be correct if the depth histogram were
     * untouched).
     *
     * Depth is capped by flattening the distribution when it exceeds
     * HUFF_MAX: halve every nonzero frequency and rebuild.  The previous
     * approach (a hand-rolled JPEG K.3 rebalance) did `cnt[i] -= 2;
     * cnt[i+1] += 1`, which loses one symbol per iteration instead of
     * preserving the count — for a bell-shaped byte distribution with 30
     * over-long codes that dropped 18 of 85 symbols, and the leftovers
     * were then all forced to 12 bits: 5.42 bits/symbol against an
     * entropy of 3.63 (48% worse than the 3.67 the tree actually gave). */
    for (int attempt = 0; attempt < 40; attempt++) {
        HPair heap[512]; int n=0;
        for(int i=0;i<256;i++) if(f[i]) hp_push(heap,&n,f[i],i);
        if(n==0) return;
        if(n==1){lens[heap[0].sym]=1;return;}
        int left[512],right[512],next_id=256;
        memset(left,-1,sizeof(left)); memset(right,-1,sizeof(right));
        while(n>1){HPair a=hp_pop(heap,&n),b=hp_pop(heap,&n);left[next_id]=a.sym;right[next_id]=b.sym;hp_push(heap,&n,a.f+b.f,next_id++);}
        int root=heap[0].sym,q[512],qh=0,qt=0,depth[512];
        memset(depth,0,sizeof(depth));
        q[qt++]=root; depth[root]=0;
        /* Leaf depths are set here from their parent's depth; re-visiting
         * a leaf in the BFS would recompute d = depth[v]+1 and add 1 to
         * every symbol's code length (a long-standing -1 bit/symbol
         * penalty). */
        while(qh<qt){int v=q[qh++];
            if(v>=256){int d=depth[v]+1;
                if(left[v]>=0){depth[left[v]]=d;q[qt++]=left[v];}
                if(right[v]>=0){depth[right[v]]=d;q[qt++]=right[v];}
            }
        }
        int maxd = 0;
        for (int i = 0; i < 256; i++) if (depth[i] > maxd) maxd = depth[i];
        if (maxd <= HUFF_MAX) {
            for (int i = 0; i < 256; i++) lens[i] = (uint8_t)depth[i];
            return;
        }
        for (int i = 0; i < 256; i++) if (f[i]) f[i] = (f[i] + 1) >> 1;
    }
    /* unreachable in practice (uniform frequencies give depth <= 8) */
    for (int i = 0; i < 256; i++) if (fq[i]) lens[i] = HUFF_MAX;
}
static void hcodes(const uint8_t lens[256],uint16_t code[256],uint8_t bits[256]){
    /* On the decode side `lens` comes straight off the wire (h_dec_p copies
     * 256 bytes of the stream into it), so it is untrusted: lengths above
     * HUFF_MAX are not a code this decoder can build, and indexing cnt[] and
     * next[], both 16 entries, with them walked off the stack.  Treat those
     * symbols as absent instead — the decoder then fails to match them and
     * returns -1, which is the answer a corrupt block deserves.  Encoder
     * callers are unaffected: hlens() never emits more than HUFF_MAX. */
    uint8_t l2[256];
    for(int i=0;i<256;i++) l2[i] = (lens[i] <= HUFF_MAX) ? lens[i] : 0;
    int cnt[16]={0};for(int i=0;i<256;i++)if(l2[i])cnt[l2[i]]++;
    uint32_t cur=0;uint16_t next[16];
    for(int L=1;L<=15;L++){next[L]=(uint16_t)cur;cur=(cur+(uint32_t)cnt[L])<<1;if(cur>65535)cur=65535;}
    for(int i=0;i<256;i++){if(l2[i]){code[i]=next[l2[i]]++;bits[i]=l2[i];}else{code[i]=0;bits[i]=0;}}
}
static uint32_t h_enc(const uint8_t *d,size_t n,const uint8_t lens[256],const uint16_t code[256],const uint8_t bits[256],uint8_t *ob,size_t cap){
    if(cap<264+n)return 0;
    for(int i=0;i<256;i++)ob[i]=lens[i];
    ob[256]=n&0xFF;ob[257]=(n>>8)&0xFF;ob[258]=(n>>16)&0xFF;ob[259]=n>>24;
    uint8_t *bst=ob+264,*op=bst;uint64_t buf=0;int cnt=0;
    /* 3-symbol unroll: max 3*15+7 = 52 bits, fits the 64-bit buffer */
    size_t i=0;
    for(;i+3<=n;i+=3){
        int s0=d[i],s1=d[i+1],s2=d[i+2];
        buf=(buf<<bits[s0])|code[s0];
        buf=(buf<<bits[s1])|code[s1];
        buf=(buf<<bits[s2])|code[s2];
        cnt+=bits[s0]+bits[s1]+bits[s2];
        while(cnt>=8){cnt-=8;*op++=(uint8_t)(buf>>cnt);buf&=((1ULL<<cnt)-1);}
    }
    for(;i<n;i++){int s=d[i];buf=(buf<<bits[s])|code[s];cnt+=bits[s];
        while(cnt>=8){cnt-=8;*op++=(uint8_t)(buf>>cnt);buf&=((1ULL<<cnt)-1);}}
    if(cnt>0)*op++=(uint8_t)(buf<<(8-cnt));
    uint32_t bsl=(uint32_t)(op-bst);ob[260]=bsl&0xFF;ob[261]=(bsl>>8)&0xFF;ob[262]=(bsl>>16)&0xFF;ob[263]=bsl>>24;
    return (uint32_t)(op-ob);
}

static int h_dec(const uint8_t **sp,const uint8_t *end,uint8_t *dst,size_t n);
static _Thread_local int t_in_worker = 0;
/* ── segment-level parallelism ─────────────────────────────────
 * v5 segments are independent: compress/decompress them across
 * threads.  Opt-in: default auto (only for enough segments),
 * CPENC_THREADS=N forces N, CPENC_NO_THREADS=1 forces 1 (PRoot's
 * scheduler can wedge under load — the env fallback stays). */
static int seg_threads(size_t n_seg) {
    /* default single-threaded: PRoot's scheduler penalizes threads,
     * and small blocks lose to thread setup.  CPENC_THREADS=N opts in
     * (real hardware / servers); CPENC_NO_THREADS=1 forces off. */
    const char *ts = getenv("CPENC_THREADS");
    if (t_in_worker) return 1;      /* never spawn threads inside a worker */
    if (getenv("CPENC_NO_THREADS")) return 1;
    if (!ts) return 1;
    int n = atoi(ts);
    if (n < 1) return 1;
    if (n > 16) n = 16;
    if ((size_t)n > n_seg) n = (int)n_seg;
    return n;
}
typedef struct { const uint8_t *src; size_t len; uint8_t *out; size_t cap; size_t outlen; int rc; int mo; } SegEncJob;
int cpenc_compress_impl(const void *sv, size_t sl, void *dv, size_t *dl,
                        int mo, int no_bwt, int no_lz);
/* Persistent workers with a shared counter.  Creating and joining a
 * thread per batch of segments cost 59-127us per thread on the phone
 * (measured), and an 8MB block has 128 segments — 8-16ms of thread
 * management against ~44ms of actual compression, which is why 3 cores
 * bought only 19%.  Now `nthr` threads are created once per block. */
typedef struct {
    const uint8_t *src; size_t sl; size_t n_seg; size_t seg_sz;
    SegEncJob *jobs; size_t next;
    pthread_mutex_t mtx; int mo;
} SegPool;
static void *seg_pool_worker(void *arg) {
    SegPool *p = arg;
    for (;;) {
        pthread_mutex_lock(&p->mtx);
        size_t k = p->next < p->n_seg ? p->next++ : (size_t)-1;
        pthread_mutex_unlock(&p->mtx);
        if (k == (size_t)-1) break;
        SegEncJob *j = &p->jobs[k];
        size_t off = k * p->seg_sz;
        size_t len = p->sl - off;
        if (len > p->seg_sz) len = p->seg_sz;
        size_t dl = j->cap;
        j->src = p->src + off; j->len = len; j->mo = p->mo;
        t_in_worker++;                 /* no nested pools */
        j->rc = cpenc_compress_impl(j->src, j->len, j->out, &dl, j->mo, 0, 1);
        t_in_worker--;
        j->outlen = dl;
    }
    return NULL;
}
typedef struct { const uint8_t *src; size_t clen; uint8_t *out; size_t cap; size_t outlen; int rc; } SegDecJob;
/* Output capacity for one segment of a segmented block.  The header's
 * original_size describes what the sub-block decodes to, EXCEPT for a
 * hex-wrapped one: there the declared size is the folded hex buffer, and the
 * text it restores to is up to twice as long.  Sizing the buffer from the
 * declaration alone made every hex segment inside a segmented block
 * undecodable - the restore did not fit and the length check then failed. */
static size_t seg_out_cap(const uint8_t *src, size_t declared) {
    if (src[3] == CPENC_STRIDE_HEX) return 2 * declared + 16;
    return declared;
}

static void *seg_dec_worker(void *arg) {
    SegDecJob *j = arg;
    size_t dl = j->cap;
    t_in_worker++;                    /* a nested segmented block must not
                                       * spawn its own pool inside a worker */
    j->rc = cpenc_decompress(j->src, j->clen, j->out, &dl);
    t_in_worker--;
    j->outlen = dl;
    return NULL;
}
/* ── const+exceptions sub-block (subtag 4): [len u32][cval u8][nex u32]
 *    [pos LEB128][val u8] each — for streams dominated by one symbol ── */
static int sub5_encode(uint8_t **po, const uint8_t *edd, const uint8_t *sym,
                       size_t n, const uint32_t cnt[256]) {
    uint8_t *out = *po;
    uint32_t dom = 0, domv = 0;
    for (int i = 0; i < 256; i++) if (cnt[i] > dom) { dom = cnt[i]; domv = (uint32_t)i; }
    if (out + 9 > edd) return -1;
    *out++ = (uint8_t)(n & 0xFF); *out++ = (uint8_t)((n >> 8) & 0xFF);
    *out++ = (uint8_t)((n >> 16) & 0xFF); *out++ = (uint8_t)((n >> 24) & 0xFF);
    *out++ = (uint8_t)domv;
    uint32_t nex = (uint32_t)n - dom;
    *out++ = (uint8_t)(nex & 0xFF); *out++ = (uint8_t)((nex >> 8) & 0xFF);
    *out++ = (uint8_t)((nex >> 16) & 0xFF); *out++ = (uint8_t)((nex >> 24) & 0xFF);
    for (size_t k = 0; k < n; k++) {
        if (sym[k] != (uint8_t)domv) {
            uint32_t pos = (uint32_t)k;
            uint8_t lb[4]; int ln = 0;
            do { lb[ln++] = (uint8_t)(pos & 0x7F); pos >>= 7;
                 if (pos) lb[ln-1] |= 0x80; } while (pos);
            for (int b2 = 0; b2 < ln; b2++) { if (out >= edd) return -1; *out++ = lb[b2]; }
            if (out >= edd) return -1;
            *out++ = sym[k];
        }
    }
    *po = out;
    return 0;
}
static int sub5_decode(const uint8_t **pp, const uint8_t *edd, uint8_t *out,
                       size_t outcap, uint32_t *plen) {
    const uint8_t *p = *pp;
    if (REMAIN(p,edd) < 9) return -1;
    uint32_t len = (uint32_t)p[0] | ((uint32_t)p[1] << 8) | ((uint32_t)p[2] << 16) | ((uint32_t)p[3] << 24); p += 4;
    if (len > outcap) return -1;
    uint8_t cval = *p++;
    uint32_t nex = (uint32_t)p[0] | ((uint32_t)p[1] << 8) | ((uint32_t)p[2] << 16) | ((uint32_t)p[3] << 24); p += 4;
    if (nex > len) return -1;
    memset(out, cval, len);
    for (uint32_t i = 0; i < nex; i++) {
        uint32_t pos = 0; int sh2 = 0;
        for (;;) {
            if (p >= edd) return -1;
            uint8_t b = *p++;
            pos |= (uint32_t)(b & 0x7F) << sh2;
            if (!(b & 0x80)) break;
            sh2 += 7;
            if (sh2 >= 28) return -1;
        }
        if (pos >= len) return -1;
        if (p >= edd) return -1;
        out[pos] = *p++;
    }
    *pp = p;
    *plen = len;
    return 0;
}
/* ── reference layer (stride 0xFB) ───────────────────────────────
 * Emit "N bytes equal what is D back" instead of the literals.  This is
 * the piece the delta/RLE tag family cannot express — RLE handles runs of
 * one value, tag5 a constant base, but a repeated *line* needs a
 * reference to an earlier region.  Logs are exactly that: a fixed
 * template repeated per record, which is why xz beats us 40x there while
 * we beat it 13x on structured tables.
 *
 * There is no dictionary: the reference points into the output produced
 * so far (head[]/chain[] are only a position index for finding
 * candidates).  The token stream is then compressed by the normal
 * machinery, so literals and tags still go through the same
 * delta/transpose/Huffman/RLE pipeline:
 *   [0][len LEB128][len bytes]                        literal run
 *   [1][len LEB128][dist LEB128]                      exact reference
 *   [2][len][dist][nex][(gap)(val)]*nex               base + patches
 * Greedy parse over a 64K-entry hash chain; distances therefore cap at
 * the chain window, which is ample for line- and record-level repeats. */
#define REF_LAZY 1
#define REF_MINM   4
/* Approximate match ("base + exceptions"): a base that is nearly, not
 * exactly, equal — the log/record case where a long body repeats with a
 * few fields changed.  Only taken when the base is long and the
 * differences sparse; on short lines (where exact matching already
 * decomposes the template into matches plus literals) merging costs more
 * than it saves and steals the anchors later matches need.  Measured:
 * 1KB records +23%, 2KB records +65%; logs and text unchanged. */
#define REF_MAXMIS  8
#define REF_APPMIN  128
#define REF_APPDIV  64
#define REF_HBITS  16
#define REF_HSZ    (1u << REF_HBITS)
#define REF_HMASK  (REF_HSZ - 1)
#define REF_MAXTRY 16   /* real data is insensitive to depth (logs, records,
                        * source all identical at 8..64); only the synthetic
                        * word-soup loses ratio, 6.49x -> 5.79x, while
                        * throughput gains 40%+ */
static uint32_t ref_h4(const uint8_t *p) {
    uint32_t v = (uint32_t)p[0] | ((uint32_t)p[1] << 8) |
                 ((uint32_t)p[2] << 16) | ((uint32_t)p[3] << 24);
    return (v * 2654435761u) >> (32 - REF_HBITS);
}
static int leb_put(uint8_t *p, uint32_t v) {
    int n = 0;
    do { uint8_t b = (uint8_t)(v & 0x7F); v >>= 7; if (v) b |= 0x80; p[n++] = b; } while (v);
    return n;
}
static int leb_len(uint32_t v) { int n = 1; while (v >= 128) { v >>= 7; n++; } return n; }
/* Exact byte size of the subtag-4 (const+exceptions) encoding of one
 * sub-stream: [len u32][cval u8][nex u32] + per exception (pos LEB)(val). */
static size_t sub5_abs_size(const uint32_t cnt[256], const uint8_t *sym, size_t n) {
    uint32_t dom = 0, domv = 0;
    for (int i = 0; i < 256; i++) if (cnt[i] > dom) { dom = cnt[i]; domv = (uint32_t)i; }
    size_t sz = 9;
    for (size_t k = 0; k < n; k++)
        if (sym[k] != (uint8_t)domv) sz += (size_t)leb_len((uint32_t)k) + 1;
    return sz;
}
/* Choose the cheapest subtag (0 const, 1 LUT, 2 Huffman, 3 RAW, 4 const+
 * exceptions) for the RLE runs/vals sub-streams and return its EXACT
 * encoded size, including the subtag byte.  The estimate and the emit
 * both call this, so a group's tag-4 bid is its real byte count — the old
 * inline estimates used floor-plus-slack and an entropy guess, and could
 * under- or overbid by thousands of bytes on long streams. */
static int sub_choose(const uint32_t cnt[256], const uint8_t *sym, size_t n, size_t *psz) {
    int nz = 0; for (int i = 0; i < 256; i++) if (cnt[i]) nz++;
    if (nz <= 1) { *psz = 7; return 0; }   /* st+len u32+val+nz-1 */
    int bits = 1; while ((1 << bits) < nz && bits < HUFF_MAX) bits++;
    size_t sz_lut = 6 + (size_t)nz + ((size_t)n * (size_t)bits + 7) / 8;
    uint8_t lens[256]; hlens(cnt, lens);
    size_t hbits = 0;
    for (int i = 0; i < 256; i++) if (cnt[i]) hbits += (size_t)cnt[i] * lens[i];
    size_t sz_h = 1 + 264 + (hbits + 7) / 8;
    size_t sz_raw = 5 + n;
    size_t sz_s5 = 1 + sub5_abs_size(cnt, sym, n);
    int st = 1; size_t best = sz_lut;
    if (sz_h < best) { best = sz_h; st = 2; }
    if (sz_raw < best) { best = sz_raw; st = 3; }
    if (sz_s5 < best) { best = sz_s5; st = 4; }
    *psz = best;
    return st;
}
/* Build the tag-4 (RLE) streams for a zero-heavy group: run lengths,
 * nonzero values, escaped 8-bit run stream, exact payload size, and the
 * chosen subtag for each sub-stream.  Selection and emit both call this,
 * so a group's RLE bid is the size it will actually occupy. */
static int rle_build(const uint8_t *gp, size_t gs,
                     uint32_t **pruns, uint8_t **pvls, size_t *pnv,
                     uint8_t **pr8, size_t *pnr8, size_t *psz,
                     int *pst_r, int *pst_v) {
    /* Sized for the worst case, not for what the caller's zero-density
     * gate happens to guarantee: the run count is bounded by the number
     * of non-zero bytes plus one, and relying on that gate made the
     * buffer overflow the moment anyone loosened it (found by an audit
     * that relaxed fq[0] > gs/2 to gs/8 and crashed here on records). */
    uint32_t *runs = malloc((gs+2)*sizeof(uint32_t));
    uint8_t *vls = malloc(gs);
    if (!runs || !vls) { free(runs); free(vls); return -1; }
    size_t nr = 0, nv = 0, i;
    if (gp[0] != 0) runs[nr++] = 0;
    for (i = 0; i < gs;) {
        if (gp[i] == 0) {
            size_t j = i; while (j < gs && gp[j] == 0) j++;
            runs[nr++] = (uint32_t)(j - i); i = j;
        } else {
            if (i > 0 && gp[i-1] != 0) runs[nr++] = 0;
            vls[nv++] = gp[i]; i++;
        }
    }
    uint8_t *r8 = malloc(nr*5+2);
    if (!r8) { free(runs); free(vls); return -1; }
    size_t nr8 = 0;
    for (size_t k = 0; k < nr; k++) {
        if (runs[k] < 255) r8[nr8++] = (uint8_t)runs[k];
        else {
            r8[nr8++] = 255;
            uint32_t r = runs[k] - 255;
            do { uint8_t b = (uint8_t)(r & 0x7F); r >>= 7;
                 if (r) b |= 0x80;
                 r8[nr8++] = b; } while (r);
        }
    }
    uint32_t rq[256] = {0}; for (size_t k = 0; k < nr8; k++) rq[r8[k]]++;
    uint32_t vq[256] = {0}; for (size_t k = 0; k < nv; k++) vq[vls[k]]++;
    size_t s_r, s_v;
    *pst_r = sub_choose(rq, r8, nr8, &s_r);
    *pst_v = sub_choose(vq, vls, nv, &s_v);
    *psz = s_r + s_v;
    *pruns = runs; *pvls = vls; *pr8 = r8;
    *pnv = nv; *pnr8 = nr8;
    return 0;
}
/* Tag 6 (bit planes): split the group into its 8 bit planes.  Each plane
 * is stored one of three ways, and the marker byte says which:
 *   0 constant — one value byte;
 *   1 packed   — pn bytes, bit i of byte i>>3 is plane bit i;
 *   2 sparse   — 4-byte count, then (LEB gap, byte) per nonzero packed
 *                byte; stores the plane's minority when it is the 1s;
 *   3 sparse, inverted — like 2, the minority is the 0s; the decoder
 *                inverts the plane after unpacking.
 * Exact payload size, filling mk[8] so the emit can reuse the decision.
 * Constant and near-constant planes are the win: quantisation scale
 * factors (E8M0) have 7 planes that barely move, while a byte Huffman
 * code must still pay for the bits that do. */
static void plane_pack(const uint8_t *gp, size_t gs, int b, uint8_t *pk) {
    size_t i = 0;
    for (; i + 8 <= gs; i += 8) {
        pk[i >> 3] = (uint8_t)(((gp[i  ] >> b) & 1)       | ((gp[i+1] >> b) & 1) << 1 |
                               ((gp[i+2] >> b) & 1) << 2  | ((gp[i+3] >> b) & 1) << 3 |
                               ((gp[i+4] >> b) & 1) << 4  | ((gp[i+5] >> b) & 1) << 5 |
                               ((gp[i+6] >> b) & 1) << 6  | ((gp[i+7] >> b) & 1) << 7);
    }
    if (i < gs) {
        uint8_t v = 0;
        for (int j = 0; i + (size_t)j < gs; j++)
            v |= (uint8_t)(((gp[i + (size_t)j] >> b) & 1) << j);
        pk[i >> 3] = v;
    }
}
static size_t bitplane_plan(const uint8_t *gp, size_t gs, const size_t bcnt[8],
                            uint8_t mk[8]) {
    size_t pn = (gs + 7) / 8;
    uint8_t *pk = malloc(pn ? pn : 1);
    if (!pk) return (size_t)-1;
    size_t sz = 8;                        /* marker bytes */
    for (int b = 0; b < 8; b++) {
        size_t o = bcnt[b];
        if (o == 0 || o == gs) { mk[b] = 0; sz += 1; continue; }
        /* A set byte costs at least 2 bytes and holds at most 8 minority
         * bits, so a plane with a dense minority cannot go sparse — decide
         * that from the counts and skip packing it (random data lands here
         * for all 8 planes). */
        {
            size_t mn0 = o < gs - o ? o : gs - o;
            if (4 + 2 * ((mn0 + 7) / 8) >= pn) { mk[b] = 1; sz += pn; continue; }
        }
        plane_pack(gp, gs, b, pk);
        int inv = (o > gs - o);           /* minority are zeros */
        if (inv) for (size_t i = 0; i < pn; i++) pk[i] ^= 0xFF;
        size_t nz = 0;
        for (size_t i = 0; i < pn; i++) nz += (pk[i] != 0);
        size_t sp = pn;
        /* every set byte costs at least 2 bytes, so below this the gaps
         * cannot possibly beat the packed form — skip the walk when they
         * cannot (packed planes are the common case) */
        if (4 + 2 * nz < pn) {
            sp = 4; uint32_t prev = 0;
            for (size_t i = 0; i < pn; i++) if (pk[i]) {
                sp += (size_t)leb_len((uint32_t)i - prev) + 1;
                prev = (uint32_t)i;
            }
        }
        if (sp < pn) { mk[b] = (uint8_t)(inv ? 3 : 2); sz += sp; }
        else { mk[b] = 1; sz += pn; }
    }
    free(pk);
    return sz;
}
/* Returns token-stream length, or 0 if the block is not worth matching. */
/* Read-only longest-match search at position i (no hash insertion), so
 * the lazy step below can look one byte ahead without disturbing the
 * chain.  Same walk as the greedy loop it replaces. */
static void ref_find(const uint8_t *s, size_t n, size_t i,
                     const int32_t *head, const int32_t *chain,
                     int *bestp, size_t *bestdp) {
    int best = 0; size_t bestd = 0;
    if (i + REF_MINM <= n) {
        int32_t c = head[ref_h4(s + i)]; int tries = 0;
        while (c >= 0 && tries < REF_MAXTRY) {
            size_t cd = i - (size_t)c;
            if (cd && s[c] == s[i] && s[c+1] == s[i+1] &&
                s[c+2] == s[i+2] && s[c+3] == s[i+3]) {
                if (best && i + best < n && s[c + best] != s[i + best]) {
                    c = chain[c & REF_HMASK]; tries++; continue;
                }
                size_t L = 0;
                while (i + L < n && s[c + L] == s[i + L]) L++;
                if ((int)L > best) { best = (int)L; bestd = cd; }
            }
            c = chain[c & REF_HMASK]; tries++;
        }
    }
    *bestp = best; *bestdp = bestd;
}
static size_t ref_encode(const uint8_t *s, size_t n, uint8_t *t, size_t cap) {
    if (n < 256) return 0;
    int32_t *head = malloc(REF_HSZ * sizeof(int32_t));
    int32_t *chain = malloc(REF_HSZ * sizeof(int32_t));
    if (!head || !chain) { free(head); free(chain); return 0; }
    for (uint32_t k = 0; k < REF_HSZ; k++) head[k] = -1;
    size_t ti = 0, i = 0, lit0 = 0;
    int fail = 0;
    while (i < n && !fail) {
        int best = 0; size_t bestd = 0;
        ref_find(s, n, i, head, chain, &best, &bestd);
        if (REF_LAZY && best >= REF_MINM && i + 1 < n) {
            /* Lazy step: if a strictly longer match starts at i+1, a
             * literal here and that match there is cheaper overall. */
            int best2 = 0; size_t bestd2 = 0;
            ref_find(s, n, i + 1, head, chain, &best2, &bestd2);
            if (best2 > best) {
                if (i + REF_MINM <= n) {
                    uint32_t hh = ref_h4(s + i);
                    chain[i & REF_HMASK] = head[hh];
                    head[hh] = (int32_t)i;
                }
                i++;
                continue;
            }
        }
        if (best >= REF_MINM) {
            /* extend the base with a few mismatches -> base + patches */
            size_t alen = 0; int ane = 0;
            uint32_t aeoff[REF_MAXMIS]; uint8_t aeval[REF_MAXMIS];
            if (bestd) {
                size_t c = i - bestd, k = 0, lastgood = 0; int ne = 0;
                while (i + k < n && c + k < i) {
                    if (s[i + k] == s[c + k]) { k++; lastgood = k; }
                    else {
                        if (ne >= REF_MAXMIS) break;
                        aeoff[ne] = (uint32_t)k; aeval[ne] = s[i + k]; ne++; k++;
                    }
                }
                while (ne > 0 && aeoff[ne - 1] >= lastgood) ne--;   /* trailing
                    * mismatches would be literals anyway — drop them */
                if (lastgood >= REF_APPMIN && (size_t)ne * REF_APPDIV < lastgood) {
                    alen = lastgood; ane = ne;
                }
            }
            long cex = 1 + leb_len((uint32_t)bestd) + leb_len((uint32_t)best);
            long capx = 0;
            if (alen) {
                capx = 1 + leb_len((uint32_t)bestd) + leb_len((uint32_t)alen) + leb_len((uint32_t)ane);
                uint32_t prev = 0;
                for (int k = 0; k < ane; k++) { capx += leb_len(aeoff[k] - prev) + 1; prev = aeoff[k]; }
            }
            int use_ap = (alen && (long)alen - capx > (long)best - cex);
            size_t mlen = use_ap ? alen : (size_t)best;

            size_t ll = i - lit0;
            if (ll) {
                if (ti + 1 + 5 + ll > cap) { fail = 1; break; }
                t[ti++] = 0;
                ti += (size_t)leb_put(t + ti, (uint32_t)ll);
                memcpy(t + ti, s + lit0, ll); ti += ll;
            }
            if (ti + 1 + 5 + 5 + 5 + (size_t)ane * 6 > cap) { fail = 1; break; }
            if (use_ap) {
                t[ti++] = 2;
                ti += (size_t)leb_put(t + ti, (uint32_t)alen);
                ti += (size_t)leb_put(t + ti, (uint32_t)bestd);
                ti += (size_t)leb_put(t + ti, (uint32_t)ane);
                uint32_t prev = 0;
                for (int k = 0; k < ane; k++) {
                    ti += (size_t)leb_put(t + ti, aeoff[k] - prev);
                    t[ti++] = aeval[k];
                    prev = aeoff[k];
                }
            } else {
                t[ti++] = 1;
                ti += (size_t)leb_put(t + ti, (uint32_t)best);
                ti += (size_t)leb_put(t + ti, (uint32_t)bestd);
            }
            for (size_t k = 0; k < mlen; k++) {
                size_t p = i + k;
                if (p + REF_MINM <= n) {
                    uint32_t hh = ref_h4(s + p);
                    chain[p & REF_HMASK] = head[hh];
                    head[hh] = (int32_t)p;
                }
            }
            i += mlen; lit0 = i;
        } else {
            if (i + REF_MINM <= n) {
                uint32_t hh = ref_h4(s + i);
                chain[i & REF_HMASK] = head[hh];
                head[hh] = (int32_t)i;
            }
            i++;
        }
    }
    if (!fail) {
        size_t ll = n - lit0;
        if (ll) {
            if (ti + 1 + 5 + ll > cap) fail = 1;
            else {
                t[ti++] = 0;
                ti += (size_t)leb_put(t + ti, (uint32_t)ll);
                memcpy(t + ti, s + lit0, ll); ti += ll;
            }
        }
    }
    free(head); free(chain);
    return fail ? 0 : ti;
}
static int ref_decode(const uint8_t *t, size_t tn, uint8_t *out, size_t outn) {
    size_t ti = 0, o = 0;
    while (ti < tn) {
        uint8_t tag = t[ti++];
        uint32_t len = 0; int sh = 0;
        for (;;) {
            if (ti >= tn || sh > 28) return -1;
            uint8_t b = t[ti++];
            len |= (uint32_t)(b & 0x7F) << sh;
            if (!(b & 0x80)) break;
            sh += 7;
        }
        if (tag == 0) {
            if (ti + len > tn || o + len > outn) return -1;
            memcpy(out + o, t + ti, len); o += len; ti += len;
        } else if (tag == 1 || tag == 2) {
            uint32_t dist = 0; sh = 0;
            for (;;) {
                if (ti >= tn || sh > 28) return -1;
                uint8_t b = t[ti++];
                dist |= (uint32_t)(b & 0x7F) << sh;
                if (!(b & 0x80)) break;
                sh += 7;
            }
            if (!dist || dist > o || o + len > outn) return -1;
            for (uint32_t k = 0; k < len; k++) out[o + k] = out[o + k - dist];
            if (tag == 2) {
                uint32_t nex = 0; sh = 0;
                for (;;) {
                    if (ti >= tn || sh > 28) return -1;
                    uint8_t b = t[ti++];
                    nex |= (uint32_t)(b & 0x7F) << sh;
                    if (!(b & 0x80)) break;
                    sh += 7;
                }
                uint32_t prev = 0;
                for (uint32_t e = 0; e < nex; e++) {
                    uint32_t gap = 0; sh = 0;
                    for (;;) {
                        if (ti >= tn || sh > 28) return -1;
                        uint8_t b = t[ti++];
                        gap |= (uint32_t)(b & 0x7F) << sh;
                        if (!(b & 0x80)) break;
                        sh += 7;
                    }
                    prev += gap;
                    if (prev >= len || ti >= tn) return -1;
                    out[o + prev] = t[ti++];
                }
            }
            o += len;
        } else return -1;
    }
    return (o == outn) ? 0 : -1;
}

int cpenc_compress_impl(const void *sv, size_t sl, void *dv, size_t *dl,
                        int mo, int no_bwt, int no_lz);


/* ─── reference layer v2: columnar sub-streams ──────────────────
 * The v1 token stream interleaves [tag][len][dist][...] and is then
 * handed to the generic compressor as ONE block, so several very
 * different statistics — the near-constant tag column, the small
 * match lengths, the wide distances, the literal text — all share a
 * single stride/order/Huffman choice.  Regrouping them into
 * homogeneous columns and compressing each with its own model is a
 * large win on source-like text: measured 2.27x -> 3.5x on 152KB of
 * C.  Distances split into a fixed u16 low half plus a varint high
 * half: the near-uniform 16-bit values let the compressor's stride
 * detection model the two byte planes separately, while a mixed-length
 * varint stream defeats that (23021 -> 21520 bytes on the corpus); the
 * high half is almost always 0, so it costs a few hundred bytes.
 * The matcher (ref_encode) is untouched; only the serialization
 * changes, and the decoder rebuilds the v1 stream so ref_decode stays
 * as it is.  v1 streams (stride CPENC_STRIDE_REF) still decode.
 * Layout: [u32 tn][u32 complen[6]][u32 uncomplen[6]][6 x cpenc block];
 * columns: tags, lens, dist-lo, dist-hi, exceptions, literals. */
#define REF2_NSUB 6
#define REF2_HDR  (4 + 8 * REF2_NSUB)
/* Split the v1 token stream into 6 columns laid out contiguously in
 * col[] (or, if col is NULL, only compute the per-column lengths).
 * Returns the total columnar size. */
static size_t ref_split2(const uint8_t *tok, size_t tn, uint8_t *col,
                         size_t plen[REF2_NSUB]) {
    size_t nt = 0, nl = 0, ndl = 0, ndh = 0, ne = 0, nx = 0, ti = 0;
    while (ti < tn) {                                  /* pass 1: sizes */
        uint8_t tag = tok[ti++]; nt++;
        uint32_t len = 0; int sh = 0;
        while (1) {
            uint8_t b = tok[ti++]; nl++;
            len |= (uint32_t)(b & 0x7F) << sh;
            if (!(b & 0x80)) break;
            sh += 7;
        }
        if (tag == 0) { nx += len; ti += len; }
        else {
            uint32_t dist = 0; sh = 0;
            while (1) {
                uint8_t b = tok[ti++];
                dist |= (uint32_t)(b & 0x7F) << sh;
                if (!(b & 0x80)) break;
                sh += 7;
            }
            ndl += 2;                                  /* dist low 16 bits */
            uint32_t dh = dist >> 16;                  /* dist high bits */
            do { ndh++; dh >>= 7; } while (dh);
            if (tag == 2) {
                uint32_t nex = 0; sh = 0;
                while (1) {
                    uint8_t b = tok[ti++]; ne++;
                    nex |= (uint32_t)(b & 0x7F) << sh;
                    if (!(b & 0x80)) break;
                    sh += 7;
                }
                for (uint32_t e = 0; e < nex; e++) {
                    while (1) { uint8_t b = tok[ti++]; ne++; if (!(b & 0x80)) break; }
                    ne++; ti++;                        /* exception value byte */
                }
            }
        }
    }
    plen[0] = nt; plen[1] = nl; plen[2] = ndl; plen[3] = ndh; plen[4] = ne; plen[5] = nx;
    if (!col) return nt + nl + ndl + ndh + ne + nx;
    uint8_t *ct = col, *cl = col + nt, *cdl = col + nt + nl,
            *cdh = col + nt + nl + ndl, *ce = col + nt + nl + ndl + ndh,
            *cx = col + nt + nl + ndl + ndh + ne;
    size_t q1 = 0, q2 = 0, q3 = 0, q4 = 0, q5 = 0, q6 = 0; ti = 0;
    while (ti < tn) {                                  /* pass 2: fill */
        uint8_t tag = tok[ti++]; ct[q1++] = tag;
        uint32_t len = 0; int sh = 0;
        while (1) {
            uint8_t b = tok[ti++];
            cl[q2++] = b;
            len |= (uint32_t)(b & 0x7F) << sh;
            if (!(b & 0x80)) break;
            sh += 7;
        }
        if (tag == 0) { memcpy(cx + q6, tok + ti, len); q6 += len; ti += len; }
        else {
            uint32_t dist = 0; sh = 0;
            while (1) {
                uint8_t b = tok[ti++];
                dist |= (uint32_t)(b & 0x7F) << sh;
                if (!(b & 0x80)) break;
                sh += 7;
            }
            cdl[q3++] = (uint8_t)dist; cdl[q3++] = (uint8_t)(dist >> 8);
            uint32_t dh = dist >> 16;
            do { uint8_t b = dh & 0x7F; dh >>= 7; if (dh) b |= 0x80; cdh[q4++] = b; } while (dh);
            if (tag == 2) {
                uint32_t nex = 0; sh = 0;
                while (1) {
                    uint8_t b = tok[ti++];
                    ce[q5++] = b;
                    nex |= (uint32_t)(b & 0x7F) << sh;
                    if (!(b & 0x80)) break;
                    sh += 7;
                }
                for (uint32_t e = 0; e < nex; e++) {
                    while (1) { uint8_t b = tok[ti++]; ce[q5++] = b; if (!(b & 0x80)) break; }
                    ce[q5++] = tok[ti++];
                }
            }
        }
    }
    return nt + nl + ndl + ndh + ne + nx;
}
/* Pack tok[0..tn) into out; returns payload size or 0 on failure. */
static size_t ref_pack2(const uint8_t *tok, size_t tn, uint8_t *out, size_t cap, int mo) {
    if (tn < 16) return 0;
    size_t plen[REF2_NSUB];
    size_t colsz = ref_split2(tok, tn, NULL, plen);
    uint8_t *col = malloc(colsz + 64);
    if (!col) return 0;
    ref_split2(tok, tn, col, plen);
    size_t off[REF2_NSUB]; { size_t a = 0; for (int i = 0; i < REF2_NSUB; i++) { off[i] = a; a += plen[i]; } }
    uint8_t *blk[REF2_NSUB] = {0}; size_t blen[REF2_NSUB] = {0};
    size_t tot = REF2_HDR; int ok = 1;
    for (int i = 0; i < REF2_NSUB && ok; i++) {
        size_t bc = cpenc_max_compressed_size(plen[i]) + 64;
        uint8_t *bb = malloc(bc);
        if (!bb) { ok = 0; break; }
        size_t cl = bc;
        if (cpenc_compress_impl(col + off[i], plen[i], bb, &cl, mo, 1, 1)) {
            free(bb); ok = 0; break;
        }
        blk[i] = bb; blen[i] = cl;
        tot += blen[i];
    }
    if (ok && tot <= cap) {
        out[0]=(uint8_t)tn; out[1]=(uint8_t)(tn>>8); out[2]=(uint8_t)(tn>>16); out[3]=(uint8_t)(tn>>24);
        for (int i = 0; i < REF2_NSUB; i++) {
            uint32_t c = (uint32_t)blen[i], u = (uint32_t)plen[i];
            out[4+4*i+0]=(uint8_t)c; out[4+4*i+1]=(uint8_t)(c>>8); out[4+4*i+2]=(uint8_t)(c>>16); out[4+4*i+3]=(uint8_t)(c>>24);
            out[28+4*i+0]=(uint8_t)u; out[28+4*i+1]=(uint8_t)(u>>8); out[28+4*i+2]=(uint8_t)(u>>16); out[28+4*i+3]=(uint8_t)(u>>24);
        }
        size_t w = REF2_HDR;
        for (int i = 0; i < REF2_NSUB; i++) { memcpy(out + w, blk[i], blen[i]); w += blen[i]; }
    } else ok = 0;
    for (int i = 0; i < REF2_NSUB; i++) free(blk[i]);
    free(col);
    return ok ? tot : 0;
}
/* Rebuild the v1 token stream (exactly tn bytes) from a packed payload.
 * Returns 0 on success. */
static int ref_unpack2(const uint8_t *pl, size_t pln, uint8_t *tok, size_t tn) {
    if (pln < REF2_HDR) return -1;
    uint32_t cl[REF2_NSUB], ul[REF2_NSUB];
    for (int i = 0; i < REF2_NSUB; i++) {
        cl[i] = (uint32_t)pl[4+4*i] | ((uint32_t)pl[4+4*i+1]<<8) | ((uint32_t)pl[4+4*i+2]<<16) | ((uint32_t)pl[4+4*i+3]<<24);
        ul[i] = (uint32_t)pl[28+4*i] | ((uint32_t)pl[28+4*i+1]<<8) | ((uint32_t)pl[28+4*i+2]<<16) | ((uint32_t)pl[28+4*i+3]<<24);
    }
    size_t off = REF2_HDR;
    uint8_t *ub[REF2_NSUB] = {0};
    int rc = 0;
    for (int i = 0; i < REF2_NSUB; i++) {
        if (off + cl[i] > pln) { rc = -1; break; }
        /* ul[] is an untrusted u32 and was used as an allocation size as-is:
         * a header asking for 4 GiB is a memory-exhaustion request, not a
         * stream.  No column can be larger than a few bytes per token. */
        if (ul[i] > tn * 4 + 64) { rc = -1; break; }
        ub[i] = malloc(ul[i] ? ul[i] : 1);
        if (!ub[i]) { rc = -2; break; }
        size_t o = ul[i];
        if (cpenc_decompress(pl + off, cl[i], ub[i], &o) != 0 || o != ul[i]) { rc = -1; break; }
        off += cl[i];
    }
    if (rc == 0) {
        uint8_t *ct = ub[0], *cln = ub[1], *cdl = ub[2], *cdh = ub[3], *ce = ub[4], *cx = ub[5];
        size_t p1 = 0, p2 = 0, p3 = 0, p4 = 0, p5 = 0, p6 = 0, ti = 0;
        /* Every index below comes out of an untrusted stream: the columns are
         * decoded blocks whose contents a corrupted container controls, and
         * nothing here was checked.  One flipped bit in a REF2 block walked
         * p6 past the end of the literals column and the copy ran off the
         * heap — fuzzing under guard pages caught it.  Each read is checked
         * against its own column's length, each write against tn. */
        while (p1 < ul[0]) {
            if (ti >= tn) { rc = -1; goto out; }
            uint8_t tag = ct[p1++]; tok[ti++] = tag;
            uint32_t len = 0; int sh = 0;
            for (;;) {
                if (p2 >= ul[1] || ti >= tn || sh > 28) { rc = -1; goto out; }
                uint8_t b = cln[p2++];
                tok[ti++] = b;
                len |= (uint32_t)(b & 0x7F) << sh;
                if (!(b & 0x80)) break;
                sh += 7;
            }
            if (tag == 0) {
                if (p6 > ul[5] || (size_t)len > ul[5] - p6 ||
                    (size_t)len > tn - ti) { rc = -1; goto out; }
                memcpy(tok + ti, cx + p6, len); ti += len; p6 += len;
            } else {
                if (ul[2] < 2 || p3 > ul[2] - 2) { rc = -1; goto out; }
                uint32_t dist = (uint32_t)cdl[p3] | ((uint32_t)cdl[p3+1] << 8); p3 += 2;
                uint32_t dh = 0; sh = 0;
                for (;;) {
                    if (p4 >= ul[3] || sh > 28) { rc = -1; goto out; }
                    uint8_t b = cdh[p4++];
                    dh |= (uint32_t)(b & 0x7F) << sh;
                    if (!(b & 0x80)) break;
                    sh += 7;
                }
                dist |= dh << 16;
                do {
                    if (ti >= tn) { rc = -1; goto out; }
                    uint8_t b = dist & 0x7F; dist >>= 7; if (dist) b |= 0x80; tok[ti++] = b;
                } while (dist);
                if (tag == 2) {
                    uint32_t nex = 0; sh = 0;
                    for (;;) {
                        if (p5 >= ul[4] || ti >= tn || sh > 28) { rc = -1; goto out; }
                        uint8_t b = ce[p5++];
                        tok[ti++] = b;
                        nex |= (uint32_t)(b & 0x7F) << sh;
                        if (!(b & 0x80)) break;
                        sh += 7;
                    }
                    for (uint32_t e = 0; e < nex; e++) {
                        for (;;) {
                            if (p5 >= ul[4] || ti >= tn) { rc = -1; goto out; }
                            uint8_t b = ce[p5++]; tok[ti++] = b;
                            if (!(b & 0x80)) break;
                        }
                        if (p5 >= ul[4] || ti >= tn) { rc = -1; goto out; }
                        tok[ti++] = ce[p5++];
                    }
                }
            }
        }
        if (ti != tn) rc = -1;
    }
out:
    for (int i = 0; i < REF2_NSUB; i++) free(ub[i]);
    return rc;
}

static uint64_t ref_h16(const uint8_t *p) {
    uint64_t a, b;
    memcpy(&a, p, 8); memcpy(&b, p + 8, 8);
    a ^= a >> 33; a *= 0xff51afd7ed558ccdULL;
    b ^= b >> 29; b *= 0xc4ceb9fe1a85ec53ULL;
    return (a ^ b) >> 16;
}
/* Fraction of 16-byte windows that repeat elsewhere in the block.  The
 * sample must be chosen by CONTENT, not by position: sampling every k-th
 * byte only sees duplicates whose distance happens to be a multiple of k,
 * so a 29KB prefix repeated at a 32.6KB distance (shared-prefix KV) scored
 * 0.00% against a true duplication of 68% and the whole block lost the
 * match layer.  Keeping one position in 16 by a hash of the bytes there
 * keeps both copies of a repeated region together whatever the distance.
 *
 * The same walk also reports `bestpos`/`bestlen`: the earlier occurrence of
 * the longest match it saw (extension capped at REF_LONG) and its length.
 * Both are used by the profitability probe below.  Sampling at `bestpos`
 * rather than at the head matters for a mixed block (varied source in front
 * of a periodic table): its head is unmatchable, a 256KB head sample scored
 * the reference layer at 17.8x the main path, the layer was skipped and the
 * block lost 114x -> 62x.  A sample drawn from the most repetitive region
 * sees what the matcher will see.  `bestlen` covers the case no sample can:
 * a repeat megabytes away (4 concatenated copies of a 1MB record block
 * score bestlen 65536 but any 256KB window sees only its own copy and
 * looks like ordinary records).
 *
 * The extension is bounded twice: it stops at REF_LONG (past that the
 * "long match" test is already decided) and is only run once per
 * REF_EXT_STRIDE bytes of block, so the walks spread over the whole block
 * instead of piling up on the first hits.  Walking every hit cost libstdc++
 * 2x — it has 36% duplicate windows whose matches run tens of KB, and every
 * one was walked; a fixed count instead missed the long matches entirely on
 * a 4MB concatenated-record block, whose first 128 hits are all short
 * template matches and whose cross-copy matches only start 1MB in. */
#define REF_LONG      (1024u)   /* a match this long is "large-scale structure" */
#define REF_EXT_STRIDE (4096u)  /* walk at most one hit per this many bytes */
static double ref_dup_frac_best(const uint8_t *s, size_t n, size_t *bestpos,
                                size_t *bestlenp) {
    if (n < 4096 || n > 0x7FFFFFF0u) { *bestpos = 0; *bestlenp = 0; return 0.0; }
    int32_t *tab = malloc(65536 * sizeof(int32_t));
    if (!tab) { *bestpos = 0; *bestlenp = 0; return 1.0; }  /* OOM: permissive */
    for (int i = 0; i < 65536; i++) tab[i] = -1;
    size_t hits = 0, tot = 0, bestlen = 0, last_ext = (size_t)-1;
    *bestpos = 0;
    for (size_t i = 0; i + 16 <= n; i += 4) {           /* 4-byte grid: repeats
        * of binary data land on it, and the scan costs a quarter as much */
        uint32_t v; memcpy(&v, s + i, 4);
        if (((v * 2654435761u) >> 28) != 0) continue;   /* 1 in 16, by content */
        uint32_t h = (uint32_t)(ref_h16(s + i) & 0xFFFF);
        int32_t p = tab[h];
        if (p >= 0 && memcmp(s + (size_t)p, s + i, 16) == 0) {
            hits++;
            if (bestlen < REF_LONG && (last_ext == (size_t)-1 ||
                                       i - last_ext >= REF_EXT_STRIDE)) {
                last_ext = i;
                size_t L = 16;
                while (L < REF_LONG && i + L < n && s[(size_t)p + L] == s[i + L]) L++;
                if (L > bestlen) { bestlen = L; *bestpos = (size_t)p; }
            }
        } else tab[h] = (int32_t)i;
        tot++;
    }
    free(tab);
    *bestlenp = bestlen;
    return tot ? (double)hits / (double)tot : 0.0;
}
/* Compressed size of the reference path over the first S bytes of s: the
 * matcher followed by the normal machinery on its token stream.  0 if the
 * sample cannot be matched at all (token stream no smaller than the input). */
static size_t ref_sample_size(const uint8_t *s, size_t S, int mo) {
    size_t tcap = S + S / 64 + 64;
    uint8_t *tok = malloc(tcap);
    if (!tok) return 0;
    size_t tn = ref_encode(s, S, tok, tcap);
    size_t tot = 0;
    if (tn && tn < S) {
        size_t bcap = cpenc_max_compressed_size(tn) + 64;
        uint8_t *bb = malloc(bcap);
        if (bb) {
            size_t cl = bcap - CPENC_HEADER_SIZE;
            if (cpenc_compress_impl(tok, tn, bb + CPENC_HEADER_SIZE, &cl, mo, 1, 1) == 0)
                tot = CPENC_HEADER_SIZE + cl;
            free(bb);
        }
        size_t pcap = cpenc_max_compressed_size(tn) + 8 * REF2_NSUB + 64;
        uint8_t *pl = malloc(pcap);
        if (pl) {
            size_t pn = ref_pack2(tok, tn, pl, pcap, mo);
            if (pn && (!tot || CPENC_HEADER_SIZE + pn < tot)) tot = CPENC_HEADER_SIZE + pn;
            free(pl);
        }
    }
    free(tok);
    return tot;
}
/* Will the reference layer beat base_len?  The duplication gate above only
 * proves repeats EXIST; on a log whose fields are random the template
 * repeats but matching still loses to the main path, and the full pass was
 * pure waste (4MB random-field log: 13.3 -> 10.1 MB/s, ratio unchanged).
 * So run the real pipeline on two nested samples and fit size(S) = a + b*S,
 * which separates the fixed cost (the first occurrence is always literals)
 * from the marginal one.  Extrapolating the raw ratio instead made a
 * periodic table look unprofitable — its whole ref stream is one 4KB
 * period, so size(S)/S shrinks as S grows and a single sample overstates
 * it (ratio 1.01 vs 0.03 for the fit).
 *
 * The threshold is loose on purpose: this gate is a speed guard, and the
 * cost of being wrong is asymmetric.  A false "skip" loses a real ratio
 * (the reference layer reaches 3013x on the periodic table); a false "run"
 * costs one pass.  Measured ratios (projected/base): periodic 0.03, text
 * 0.00, log template 0.06, mixed 0.02, records 1.10, random-field log 1.91
 * — so 1.4 sits in an empty band, not against a cluster.
 *
 * A long match short-circuits the probe: it is large-scale structure no
 * sample of this size can see (it may be megabytes away), and every corpus
 * with one profited hugely — the periodic table, the single-line log, the
 * concatenated record block.  The short-match corpora are exactly the ones
 * the probe is for. */
#define REF_PROBE  (128u << 10)                /* samples: 128KB and 256KB */
static int ref_profitable(const uint8_t *s, size_t sl, size_t base_len,
                          int mo, size_t bestpos, size_t bestlen) {
    size_t S1 = REF_PROBE, S2 = REF_PROBE * 2;
    if (bestlen >= REF_LONG) return 1;
    if (S2 >= sl) return 1;                    /* small block: the pass is cheap */
    size_t start = bestpos;
    if (start + S2 > sl) start = sl - S2;
    size_t t1 = ref_sample_size(s + start, S1, mo);
    /* The sample could not be matched at all: its repetition lives further
     * away than 128KB reaches (a table whose period is 2MB), so the sample
     * is unrepresentative — defer to the duplication gate. */
    if (!t1) return 1;
    size_t t2 = ref_sample_size(s + start, S2, mo);
    if (!t2) return 1;
    double b = (t2 > t1) ? (double)(t2 - t1) / (double)(S2 - S1) : 0.0;
    double a = (double)t1 - b * (double)S1;
    if (a < 0.0) a = 0.0;
    double proj = a + b * (double)sl;
    return proj < (double)base_len * 1.4;
}
/* Try the match layer against base_len (the best size so far).  On a win
 * it writes the whole block into d and returns its length, else 0.
 * Called from two places: the RAW passthrough needs it too, because
 * records that differ in a few scattered bytes look like noise to
 * delta/BWT and fall back to RAW — returning there meant the candidate
 * never ran and 1KB-record data stored at 1.00x instead of 57x. */
static size_t ref_try(const uint8_t *s, size_t sl, uint8_t *d, size_t cap,
                     size_t base_len, int mo) {
    if (sl < 4096) return 0;
    /* No "main path is already good enough" shortcut here.  It cost a
     * periodic table 31x: the main path reached 95x and the reference
     * layer was never given a chance, when matching the 4KB period takes
     * it to 3013x.  The duplication sample below is the right gate —
     * it looks at the data, not at what the main path scored — and on
     * records, text and random data the outcome is unchanged. */
    size_t tcap = sl + sl / 64 + 64;
    uint8_t *tok = malloc(tcap);
    if (!tok) return 0;
    {
        /* Two data-driven gates, both over the whole block.  The first is a
         * spread duplication test (see ref_dup_frac_best); an earlier
         * version also encoded the first 32KB here, but that sample cannot
         * see repeats living further in — layer-0 KV is a per-token lookup
         * table whose rows repeat 2.3x and whose first 85 rows are distinct
         * — and it cost 16% of compression time on FP4 weights, where the
         * spread test says no before a byte is matched.  Matches shorter
         * than the 16-byte window are worth ~3 bytes each, so nothing of
         * value is behind that gate.  The second predicts whether matching
         * can beat the main path at all; it is the one that pays for
         * itself on data where duplication exists but matching still
         * loses. */
        if (sl < 4096) { free(tok); return 0; }
        size_t bestpos = 0, bestlen = 0;
        if (ref_dup_frac_best(s, sl, &bestpos, &bestlen) < REF_DUP_MIN) { free(tok); return 0; }
        /* Repeats exist, but can the matcher actually profit?  Predicting
         * that is cheaper than finding out: the sample probe costs a
         * fraction of a pass, and on the random-field log it saves the
         * whole pass (13.4 vs 9.9 MB/s, identical output). */
        if (!ref_profitable(s, sl, base_len, mo, bestpos, bestlen)) { free(tok); return 0; }
    }
    size_t tn = ref_encode(s, sl, tok, tcap);
    size_t res = 0;
    if (tn && tn < sl) {
        /* v1: the interleaved token stream as one block — wins on logs,
         * records and near-constant text, where matches are long and few
         * and the fixed columnar header dominates. */
        size_t bcap = cpenc_max_compressed_size(tn) + 64;
        uint8_t *bb = malloc(bcap);
        if (bb) {
            size_t cl = bcap - CPENC_HEADER_SIZE;
            if (cpenc_compress_impl(tok, tn, bb + CPENC_HEADER_SIZE, &cl, mo, 1, 1) == 0) {
                size_t tot = CPENC_HEADER_SIZE + cl;
                if (tot < base_len && tot <= cap) {
                    hdr_put((uint8_t *)bb, CPENC_STRIDE_REF, (uint32_t)sl, (uint32_t)cl);
                    memcpy(d, bb, tot);
                    res = tot;
                }
            }
            free(bb);
        }
        /* v2: columnar sub-streams — wins on source-like text, where the
         * token stream is dense short matches and each column gets its
         * own model.  Kept only if strictly smaller than v1. */
        size_t pcap = cpenc_max_compressed_size(tn) + 8 * REF2_NSUB + 64;
        uint8_t *pl = malloc(pcap);
        if (pl) {
            size_t pn = ref_pack2(tok, tn, pl, pcap, mo);
            if (pn) {
                size_t tot = CPENC_HEADER_SIZE + pn;
                size_t lim = res ? res : base_len;
                if (tot < lim && tot <= cap) {
                    hdr_put((uint8_t *)d, CPENC_STRIDE_REF2, (uint32_t)sl, (uint32_t)pn);
                    memcpy(d + CPENC_HEADER_SIZE, pl, pn);
                    res = tot;
                }
            }
            free(pl);
        }
    }
    free(tok);
    return res;
}

/* Segmented Huffman: split the symbol stream into SEGS independent
 * Huffman blocks (own table each).  Header: [SEGS][SEGS x u32 len].  Lets
 * the decoder run segments in parallel; table overhead ~SEGS*264B. */
#define HUFF_SEGS 4
static uint32_t h_enc_seg(const uint8_t *d, size_t n, uint8_t *ob, size_t cap) {
    if (cap < 1 + 4*HUFF_SEGS + 264u + n) return 0;
    size_t segn = (n + HUFF_SEGS - 1) / HUFF_SEGS;
    uint8_t *lp = ob;
    *lp++ = HUFF_SEGS;
    uint8_t *lens_out = lp;                       /* lengths, patched later */
    memset(lens_out, 0, 4*HUFF_SEGS);
    lp += 4*HUFF_SEGS;
    uint32_t total = 0;
    for (int k = 0; k < HUFF_SEGS; k++) {
        size_t off = (size_t)k * segn;
        size_t ln = (off + segn <= n) ? segn : n - off;
        if (!ln) { lens_out[k] = 0; continue; }
        uint32_t fq[256];
        freq256(d + off, ln, fq);
        uint8_t lens[256]; uint16_t code[256]; uint8_t bits[256];
        hlens(fq, lens); hcodes(lens, code, bits);
        uint32_t used = h_enc(d + off, ln, lens, code, bits, lp, cap - (size_t)(lp - ob));
        if (!used) return 0;
        lens_out[4*k+0] = (uint8_t)used;
        lens_out[4*k+1] = (uint8_t)(used >> 8);
        lens_out[4*k+2] = (uint8_t)(used >> 16);
        lens_out[4*k+3] = (uint8_t)(used >> 24);
        lp += used;
        total += used;
    }
    for (int k = 0; k < HUFF_SEGS; k++) {
        size_t off = (size_t)k * segn;
        size_t ln = (off + segn <= n) ? segn : n - off;
        (void)ln;
    }
    return (uint32_t)(lp - ob);
}
/* Exact size of the Huffman payload emitted for a whole group: h_enc, or
 * h_enc_seg once the group crosses 1MB.  Computing the size from the code
 * lengths avoids an entropy estimate whose error (even 1%) is thousands of
 * bytes on an 8MB group. */
static size_t huff_size_exact(const uint8_t *p, size_t n, const uint32_t fq[256]) {
    if (n >= (1u << 20)) {
        size_t segn = (n + HUFF_SEGS - 1) / HUFF_SEGS;
        size_t total = 1 + 4 * HUFF_SEGS;
        for (int k = 0; k < HUFF_SEGS; k++) {
            size_t off = (size_t)k * segn;
            size_t ln = (off + segn <= n) ? segn : n - off;
            if (!ln) continue;
            uint32_t f2[256]; freq256(p + off, ln, f2);
            uint8_t lens[256]; hlens(f2, lens);
            size_t hb = 0;
            for (int i = 0; i < 256; i++) if (f2[i]) hb += (size_t)f2[i] * lens[i];
            total += 264 + (hb + 7) / 8;
        }
        return total;
    }
    uint8_t lens[256]; hlens(fq, lens);
    size_t hb = 0;
    for (int i = 0; i < 256; i++) if (fq[i]) hb += (size_t)fq[i] * lens[i];
    return 264 + (hb + 7) / 8;
}
typedef struct {
    const uint8_t *data[2];  /* segment starts */
    const uint8_t *end[2];   /* segment ends */
    uint8_t *dst[2];
    size_t n[2];
    int rc;
} SegJob;
static void *seg_worker(void *arg) {
    SegJob *j = (SegJob *)arg;
    j->rc = 0;
    for (int k = 0; k < 2; k++) {
        if (!j->n[k]) continue;
        const uint8_t *seg = j->data[k];
        if (h_dec(&seg, j->end[k], j->dst[k], j->n[k])) { j->rc = -1; break; }
    }
    return NULL;
}
static int h_dec_seg(const uint8_t **sp, const uint8_t *end, uint8_t *dst, size_t n) {
    const uint8_t *src = *sp;
    if (REMAIN(src,end) < 1 + 4*HUFF_SEGS) return -1;
    int segs = *src++;
    if (segs != HUFF_SEGS) return -1;
    uint32_t lens[4];
    for (int k = 0; k < HUFF_SEGS; k++)
        lens[k] = (uint32_t)src[4*k] | ((uint32_t)src[4*k+1] << 8) |
                  ((uint32_t)src[4*k+2] << 16) | ((uint32_t)src[4*k+3] << 24);
    src += 4*HUFF_SEGS;
    size_t segn = (n + HUFF_SEGS - 1) / HUFF_SEGS;
    const uint8_t *segp[HUFF_SEGS];
    segp[0] = src;
    /* Validate every segment as it is laid out, not just the last one: the
     * lengths come from the stream, and with a raised lens[0] the later
     * starts land past `end`, where the old end-minus-start check wrapped
     * and passed.  That is how a corrupt block reached h_dec with an end
     * pointer beyond its buffer. */
    for (int k = 0; k < HUFF_SEGS; k++) {
        if (segp[k] > end || REMAIN(segp[k], end) < lens[k]) return -1;
        if (k + 1 < HUFF_SEGS) segp[k+1] = segp[k] + lens[k];
    }
    /* 2 workers x 2 segments, segments independent.  Threads are the
     * default (real hardware scales ~linearly); CPENC_NO_THREADS=1 falls
     * back to sequential where thread scheduling misbehaves (PRoot). */
    SegJob jobs[2];
    pthread_t th[2];
    for (int k = 0; k < HUFF_SEGS; k++) {
        size_t off = (size_t)k * segn;
        size_t ln = (off + segn <= n) ? segn : n - off;
        int w = k / 2, s = k % 2;
        jobs[w].data[s] = segp[k];
        jobs[w].end[s] = segp[k] + lens[k];
        jobs[w].dst[s] = dst + off;
        jobs[w].n[s] = ln;
        jobs[w].rc = 0;
    }
    if (getenv("CPENC_NO_THREADS")) {
        for (int w = 0; w < 2; w++) {
            const uint8_t *s0 = jobs[w].data[0];
            if (jobs[w].n[0] && h_dec(&s0, jobs[w].end[0], jobs[w].dst[0], jobs[w].n[0])) return -1;
            const uint8_t *s1 = jobs[w].data[1];
            if (jobs[w].n[1] && h_dec(&s1, jobs[w].end[1], jobs[w].dst[1], jobs[w].n[1])) return -1;
        }
    } else {
        if (pthread_create(&th[0], NULL, seg_worker, &jobs[0]) != 0) return -1;
        if (pthread_create(&th[1], NULL, seg_worker, &jobs[1]) != 0) return -1;
        pthread_join(th[0], NULL);
        pthread_join(th[1], NULL);
        if (jobs[0].rc || jobs[1].rc) return -1;
    }
    *sp = segp[HUFF_SEGS-1] + lens[HUFF_SEGS-1];
    return 0;
}

static int h_dec(const uint8_t **sp,const uint8_t *end,uint8_t *dst,size_t n);
static int h_dec_p(const uint8_t **sp,const uint8_t *end,uint8_t *dst,size_t n,int prefix){
    const uint8_t *src=*sp;if(REMAIN(src,end)<264)return-1;
    uint8_t lens[256];for(int i=0;i<256;i++)lens[i]=src[i];
    uint32_t ul=(uint32_t)src[256]|((uint32_t)src[257]<<8)|((uint32_t)src[258]<<16)|((uint32_t)src[259]<<24);
    uint32_t bsl=(uint32_t)src[260]|((uint32_t)src[261]<<8)|((uint32_t)src[262]<<16)|((uint32_t)src[263]<<24);
    if(bsl > (uint32_t)(end-src)-264) return -1;
    src+=264;const uint8_t *bse=src+bsl;
    uint16_t code[256];uint8_t bits[256];hcodes(lens,code,bits);
    uint16_t fast[HUFF_TABSZ];memset(fast,0,sizeof(fast));
    for(int s=0;s<256;s++){if(lens[s]&&lens[s]<=HUFF_MAX){uint16_t c=code[s];int b=lens[s];
        int xm=1<<(HUFF_MAX-b);
        uint16_t v=(uint16_t)(((uint32_t)b<<8)|s);
        for(int x=0;x<xm;x++){uint32_t fi=((uint32_t)c<<(HUFF_MAX-b))|(uint32_t)x;if(fi<HUFF_TABSZ)fast[fi]=v;}}}
    size_t tgt=ul<n?ul:n;uint64_t buf=0;int cnt=0;
    int norefill = getenv("CPENC_NOREFILL") ? 1 : 0;
    for(size_t dc=0;dc<tgt;){
        /* batch refill: one 4-byte load covers ~3-4 symbols instead of
         * per-symbol byte reads.  The bitstream is MSB-first, so each
         * loaded byte must be bit-reversed before shifting into the
         * window (a plain LE load would feed LSB-first bits to the
         * consumer).  cnt stays <= 64 (refill gated at 32). */
        if(!norefill&&cnt<=32&&src+4<=bse){
            uint32_t w;memcpy(&w,src,4);src+=4;
            /* old bits move up (like the per-byte buf<<8), new bytes land
             * in the low end — MSB-first consumer reads the old bits first */
            buf=(buf<<32)|(uint64_t)__builtin_bswap32(w);
            cnt+=32;
        }
        while(cnt<HUFF_MAX&&src<bse){buf=(buf<<8)|(*src++);cnt+=8;}
        if(cnt<HUFF_MAX){buf<<=(HUFF_MAX-cnt);cnt=HUFF_MAX;}
        uint16_t v=fast[(buf>>(cnt-HUFF_MAX))&(HUFF_WSZ-1)];
        if(!(v>>8)) return -1;
        if((v>>8)==1){
            uint8_t s1=(uint8_t)fast[0];  /* code=0 1-bit symbol (not v: v may be the code=1 one) */
            uint16_t v2=fast[HUFF_WSZ/2];
            uint8_t s2=(uint8_t)v2;
            int s2len=v2>>8;
            for(;;){
                if(cnt<=1){
                    if(src<bse){buf=(buf<<8)|(*src++);cnt+=8;}
                    else{buf<<=(HUFF_MAX-cnt);cnt=HUFF_MAX;break;}
                }
                if((buf>>(cnt-1))&1){
                    if(s2len!=1) break;
                    dst[dc]=(uint8_t)((uint16_t)s2+(dc?(uint16_t)dst[dc-1]:0));
                    dc++;cnt-=1;
                    if(dc>=tgt)break;
                    continue;
                }
                dst[dc]=(uint8_t)((uint16_t)s1+(dc?(uint16_t)dst[dc-1]:0));
                dc++;cnt-=1;
                if(dc>=tgt)break;
            }
            continue;
        }
        dst[dc]=(uint8_t)((uint16_t)v+(dc?(uint16_t)dst[dc-1]:0));
        cnt-=(v>>8);dc++;
    }
    *sp=bse;return 0;
}
static int h_dec(const uint8_t **sp,const uint8_t *end,uint8_t *dst,size_t n){
    const uint8_t *src=*sp;if(REMAIN(src,end)<264)return-1;
    uint8_t lens[256];for(int i=0;i<256;i++)lens[i]=src[i];
    uint32_t ul=(uint32_t)src[256]|((uint32_t)src[257]<<8)|((uint32_t)src[258]<<16)|((uint32_t)src[259]<<24);
    uint32_t bsl=(uint32_t)src[260]|((uint32_t)src[261]<<8)|((uint32_t)src[262]<<16)|((uint32_t)src[263]<<24);
    /* bsl comes from the stream: reject out-of-range before any pointer
     * arithmetic, or a corrupted block reads past the input buffer */
    if(bsl > (uint32_t)(end-src)-264) return -1;
    src+=264;const uint8_t *bse=src+bsl;
    uint16_t code[256];uint8_t bits[256];hcodes(lens,code,bits);
    /* fast table: high byte = code length, low byte = symbol — one lookup
     * yields both, avoiding the second lens[] access per symbol */
    uint16_t fast[HUFF_TABSZ];memset(fast,0,sizeof(fast));
    for(int s=0;s<256;s++){if(lens[s]&&lens[s]<=HUFF_MAX){uint16_t c=code[s];int b=lens[s];
        int xm=1<<(HUFF_MAX-b);
        uint16_t v=(uint16_t)(((uint32_t)b<<8)|s);
        for(int x=0;x<xm;x++){uint32_t fi=((uint32_t)c<<(HUFF_MAX-b))|(uint32_t)x;if(fi<HUFF_TABSZ)fast[fi]=v;}}}
    size_t tgt=ul<n?ul:n;uint64_t buf=0;int cnt=0;
    int norefill = getenv("CPENC_NOREFILL") ? 1 : 0;
    for(size_t dc=0;dc<tgt;){
        /* batch refill: one 4-byte load covers ~3-4 symbols instead of
         * per-symbol byte reads.  The bitstream is MSB-first, so each
         * loaded byte must be bit-reversed before shifting into the
         * window (a plain LE load would feed LSB-first bits to the
         * consumer).  cnt stays <= 64 (refill gated at 32). */
        if(!norefill&&cnt<=32&&src+4<=bse){
            uint32_t w;memcpy(&w,src,4);src+=4;
            /* old bits move up (like the per-byte buf<<8), new bytes land
             * in the low end — MSB-first consumer reads the old bits first */
            buf=(buf<<32)|(uint64_t)__builtin_bswap32(w);
            cnt+=32;
        }
        while(cnt<HUFF_MAX&&src<bse){buf=(buf<<8)|(*src++);cnt+=8;}
        if(cnt<HUFF_MAX){buf<<=(HUFF_MAX-cnt);cnt=HUFF_MAX;}
        uint16_t v=fast[(buf>>(cnt-HUFF_MAX))&(HUFF_WSZ-1)];
        if(!(v>>8)) return -1;
        /* fast path: canonical 1-bit codes (top bit 0) — skewed streams
         * like FP4 deltas are dominated by one symbol; consume runs of it
         * at ~2 ops/symbol instead of the full table lookup */
        if((v>>8)==1){
            uint8_t s1=(uint8_t)fast[0];  /* code=0 1-bit symbol (not v: v may be the code=1 one) */
            uint16_t v2=fast[HUFF_WSZ/2];
            uint8_t s2=(uint8_t)v2;
            int s2len=v2>>8;
            for(;;){
                if(cnt<=1){
                    if(src<bse){buf=(buf<<8)|(*src++);cnt+=8;}
                    else{buf<<=(HUFF_MAX-cnt);cnt=HUFF_MAX;break;}
                }
                if((buf>>(cnt-1))&1){
                    /* other symbol: 1-bit only if it really is; longer
                     * codes go to the slow path (which won't bounce back) */
                    if(s2len!=1) break;
                    dst[dc++]=s2;cnt-=1;
                    if(dc>=tgt)break;
                    continue;
                }
                dst[dc++]=s1;cnt-=1;
                if(dc>=tgt)break;
            }
            continue;
        }
        dst[dc]=(uint8_t)v;cnt-=(v>>8);
        dc++;
    }
    *sp=bse;return 0;
}

/* ── RLE2: bzip2-style zero-run coding on the MTF stream (v7) ─────
 * MTF turns BWT clustering into long runs of rank 0, which the generic
 * order-0 Huffman backend still pays a full symbol each for.  RLE2 folds a
 * run of z+1 zeros into ceil(log2(z+2)) bytes, bzip2's scheme:
 *
 *   code 0 = RUNA (run bit 0), 1 = RUNB (run bit 1),
 *   2..254 = a literal byte 1..253 (value = code - 1),
 *   255 = escape followed by the raw byte.
 *
 * The escape exists because MTF rank 255 would otherwise need code 256.
 * Ranks 254 and 255 therefore cost two bytes; every other literal costs
 * one, and runs cost far less.  Same scheme as the /tmp/algo/codec.c
 * prototype, which measured +8% ratio on text at 512K blocks.
 *
 * rle2_enc returns the encoded length, or 0 if the output does not fit the
 * capacity (0 can never be a valid encoding: the input is never empty on
 * this path).  rle2_dec returns the number of bytes it produced — compare
 * against the expected length — or 0 on a malformed stream. */
static size_t rle2_enc(const uint8_t *in, size_t n, uint8_t *out, size_t cap) {
    size_t o = 0, i = 0;
    while (i < n) {
        if (in[i] == 0) {
            size_t j = i; while (j < n && in[j] == 0) j++;
            uint32_t z = (uint32_t)(j - i) - 1;
            do {
                if (o + 2 > cap) return 0;
                out[o++] = (uint8_t)((z & 1) ? 1 : 0);
                if (z < 2) break;
                z = (z - 2) / 2;
            } while (1);
            i = j;
        } else {
            if (o + 2 > cap) return 0;
            uint8_t b = in[i++];
            if (b <= 253) out[o++] = (uint8_t)(b + 1);
            else { out[o++] = 255; out[o++] = b; }
        }
    }
    return o;
}
static size_t rle2_dec(const uint8_t *in, size_t m, uint8_t *out, size_t n) {
    size_t i = 0, k = 0;
    while (i < n) {
        if (k >= m) return 0;
        uint8_t b = in[k++];
        if (b <= 1) {
            int32_t es = -1; uint32_t N = 1; int have = 0;
            for (;;) {
                es += (b == 0) ? (int32_t)N : (int32_t)(2 * N);
                N <<= 1;
                if (k >= m) break;
                b = in[k++];
                if (b > 1) { have = 1; break; }
            }
            es++;
            if (es < 0 || (size_t)es > n - i) return 0;
            memset(out + i, 0, (size_t)es);
            i += (size_t)es;
            if (!have) { if (i != n) return 0; break; }
        }
        if (b == 255) {
            if (k >= m || i >= n) return 0;
            out[i++] = in[k++];
        } else {
            if (b < 2) return 0;
            out[i++] = (uint8_t)(b - 1);
        }
    }
    return i;
}

/* Fused RLE2 + MTF decode, for the v7 BWT2 path.  The two passes it
 * replaces touch a n1-byte intermediate buffer that nothing else reads, and
 * the run case has a shortcut: a run of zeros decodes to MTF rank 0, whose
 * move-to-front is a no-op, so the whole run is the current front byte
 * repeated (a memset) and the MTF state does not change.  Only literals pay
 * the O(rank) shift.  Saves ~1.8 ms/MB on 512K text blocks.
 * Same contract as rle2_dec: returns the number of bytes produced, or 0 on
 * a malformed stream. */
static size_t rle2_mtf_dec(const uint8_t *in, size_t m, uint8_t *out, size_t n) {
    uint8_t st[256];
    for (int i = 0; i < 256; i++) st[i] = (uint8_t)i;
    size_t i = 0, k = 0;
    while (i < n) {
        if (k >= m) return 0;
        uint8_t b = in[k++];
        if (b <= 1) {
            int32_t es = -1; uint32_t N = 1; int have = 0;
            for (;;) {
                es += (b == 0) ? (int32_t)N : (int32_t)(2 * N);
                N <<= 1;
                if (k >= m) break;
                b = in[k++];
                if (b > 1) { have = 1; break; }
            }
            es++;
            if (es < 0 || (size_t)es > n - i) return 0;
            memset(out + i, st[0], (size_t)es);   /* rank 0 is a no-op move */
            i += (size_t)es;
            if (!have) { if (i != n) return 0; break; }
        }
        uint8_t r;
        if (b == 255) {
            if (k >= m || i >= n) return 0;
            r = in[k++];
        } else {
            if (b < 2) return 0;
            r = (uint8_t)(b - 1);
        }
        uint8_t v = st[r];
        out[i++] = v;
        for (int j = r; j > 0; j--) st[j] = st[j - 1];
        st[0] = v;
    }
    return i;
}

/* ═══════════════════════════════════════════════════════════════
   Public API
   ═══════════════════════════════════════════════════════════════ */

size_t cpenc_max_compressed_size(size_t sl) {
    return CPENC_HEADER_SIZE+16+sl+sl/2+16384;
}

/* v7: the BWT2 inner call needs the main path with the lag-delta
 * preprocessor disabled.  The lag scan (v3) is chosen on a sample and is
 * right for framed records, but on an RLE2-coded MTF stream it fires often
 * and destroys the run-code structure: measured on the 16MB text corpus it
 * cost 224KB (5.76x -> 5.50x), while the zero-dense MTF stream it was
 * already tuned to skip.  `no_lag` is the only difference; every other
 * caller (including the recursive ones inside) goes through the unchanged
 * 7-argument entry point, so no existing path changes behaviour. */
static int cpenc_compress_core(const void *sv, size_t sl, void *dv, size_t *dl,
                               int mo, int no_bwt, int no_lz, int no_lag);
int cpenc_compress_impl(const void *sv, size_t sl, void *dv, size_t *dl, int mo, int no_bwt, int no_lz) {
    return cpenc_compress_core(sv, sl, dv, dl, mo, no_bwt, no_lz, 0);
}
static int cpenc_compress_core(const void *sv, size_t sl, void *dv, size_t *dl,
                               int mo, int no_bwt, int no_lz, int no_lag) {
    lg2_init();
    const uint8_t *s=(const uint8_t*)sv; uint8_t *d=(uint8_t*)dv;
    if(*dl<CPENC_HEADER_SIZE) return -1;
    if(sl>0xFFFFFFFFU) return -2;
    mo=mo<0?0:mo>4?4:mo;
    if(sl==0){memset(d,0,CPENC_HEADER_SIZE);hdr_put((uint8_t*)d,CPENC_STRIDE_RAW,0,0);
        *dl=CPENC_HEADER_SIZE;return 0;}
    /* Tiny inputs: the stride/group machinery breaks down (gs can be 0)
     * and compression gains are nil — store RAW, always safe to decode. */
    if(sl<64){
        if(*dl<CPENC_HEADER_SIZE+1+sl) return -1;
        hdr_put((uint8_t*)d,CPENC_STRIDE_RAW,(uint32_t)sl,(uint32_t)(1+sl));
        d[CPENC_HEADER_SIZE]=CPENC_LAG_NONE;
        memcpy(d+CPENC_HEADER_SIZE+1,s,sl);*dl=CPENC_HEADER_SIZE+1+sl;return 0;
    }

    /* Algebraic gate: w is the f(x1)=f(x2) statistic over the window
     * vector.  A block with almost no repeated 4-byte windows cannot win
     * with a match layer, so REF/BWT are skipped instead of trial-
     * compressed.  Measured: floats and BF16 lose zero ratio while the
     * encoder gets 1.6-2.4x faster; text/records/ELF keep w >= 0.6 and run
     * both candidates unchanged.  Inner recursions (no_bwt && no_lz) skip
     * the measurement; CPENC_NO_GATE disables the gate for A/B tests. */
    int allow_ref = 1, allow_bwt = 1;
    if (!(no_bwt && no_lz) && (!getenv("CPENC_NO_GATE") || getenv("CPENC_VFEAT"))) {
        VFeat vf; vfeat(s, sl, &vf, getenv("CPENC_VFEAT") || getenv("CPENC_CHOOSE"));
        if (getenv("CPENC_VFEAT"))
            fprintf(stderr, "VFEAT n=%zu z=%.3f a2=%.3f a4=%.3f a8=%.3f w=%.4f "
                    "r1=%.3f r2=%.3f r4=%.3f r8=%.3f r16=%.3f kbest=%d drho=%.3f lr=%.4f lr16=%.4f lr32=%.4f cover=%.3f ap1=%.3f ap2=%.3f rsim=%.3f rdup=%.3f dapprox=%.3f dshape=%.3f ms=%.3f ml=%.3f\n",
                    sl, vf.z, vf.a2, vf.a4, vf.a8, vf.w,
                    vf.rho[1], vf.rho[2], vf.rho[4], vf.rho[8], vf.rho[16],
                    vf.kbest, vf.drho, vf.lr, vf.lr16, vf.lr32, vf.cover, vf.ap1, vf.ap2, vf.rsim, vf.rdup, vf.dapprox, vf.dshape, vf.ml_short, vf.ml_long);
        allow_ref = vf.w >= 0.005;
        allow_bwt = vf.w >= 0.10;
        /* Threshold rule (fitted on 32 text blocks, 30/32): when both
         * approximation rates are high the block is REF territory, else
         * BWT.  Exclusive choice = one match layer instead of two. */
        if (getenv("CPENC_CHOOSE") && allow_ref && allow_bwt) {
            /* REF -> BWT -> entropy.  Two REF regimes measured: smooth
             * near-duplicates (both approximation rates high) and
             * fixed-width records (k-rate peak at 8, as in ELF).  Blocks
             * near the boundary run both candidates. */
            int ref_conf = (vf.ap1 >= 0.241 && vf.ap2 >= 0.129) ||
                           (vf.kbest == 8 && vf.drho >= 0.10);
            int bwt_conf = (vf.ap1 <= 0.200 || vf.ap2 <= 0.085);
            if (ref_conf) allow_bwt = 0;
            else if (bwt_conf) allow_ref = 0;
        }
    }

    /* v4: hex-ascii layer — text-form numbers (log hex dumps) waste the
     * delta machinery on half-byte interleave ('F' 'F' = 0xFF, two chars
     * per value).  Detect hex-heavy text, fold 2 chars -> 1 byte, and keep
     * the non-hex layout chars (newlines etc.) in a position table so the
     * byte stream is restored exactly.  Encoded as [bin_len u32][n_sep u32]
     * [(idx LEB128)(char)]*n_sep [bin bytes]. */
    int hex_mode = 0;
    uint8_t *hexbuf = NULL;
    /* Sampling gate first: deciding whether to try the hex layer by
     * scanning every byte costs 12-18ms on an 8MB block — 18% of total
     * compression time, and 39% on a small KV block.  Binary data has
     * 22/256 = 8.6% hex characters and hex text is >=90%, so four spread
     * windows separate them; anything in between still gets the exact
     * scan, so the only risk is a hex file whose every sampled window
     * looks binary. */
    int hex_scan = (sl < 256) ? 0 : 1;
    if (hex_scan && sl >= 8192) {
        hex_scan = 0;
        for (int k = 0; k < 4 && !hex_scan; k++) {
            size_t off = (sl - 2048) * (size_t)k / 3;
            size_t hx = 0;
            for (size_t i = 0; i < 2048; i++) {
                uint8_t c = s[off + i];
                if ((c >= '0' && c <= '9') || (c >= 'a' && c <= 'f') ||
                    (c >= 'A' && c <= 'F')) hx++;
            }
            if (hx * 2 >= 2048) hex_scan = 1;      /* >=50%: look properly */
        }
    }
    if (hex_scan) {
        size_t hexc = 0, badc = 0, hi_upper = 0, hi_lower = 0;
        for (size_t i = 0; i < sl; i++) {
            uint8_t c = s[i];
            if ((c >= '0' && c <= '9') || (c >= 'a' && c <= 'f') || (c >= 'A' && c <= 'F')) {
                hexc++;
                if (c >= 'a') hi_lower++;
                else if (c >= 'A') hi_upper++;
            }
            else if (c == '\n' || c == '\r' || c == '\t' || c == ' ' || c == ':' || c == ',' || c == '-' || c == '|') { /* layout chars ok */ }
            else badc++;
        }
        /* >=90% hex chars, <1% foreign, uniform case ('B' vs 'b' both map to
         * one nibble — mixed case would lose letters on restore) */
        if (hexc * 10 >= sl * 9 && badc * 100 < sl && (hi_upper == 0 || hi_lower == 0)) {
            uint8_t hex_case = hi_upper > 0 ? 1 : 0;
            /* A lone nibble (odd digit count, or one left before a foreign
             * char) travels in the separator table as a character, and it
             * has to carry the block's case: writing it lowercase turned
             * the last 'A' of 257 of them into 'a' on the way back. */
            uint8_t hex_lit = (uint8_t)(hex_case ? 'A' : 'a');
            /* The separator table is written at offset 13+sl and moved down
             * afterwards, so it needs its own worst-case scratch.  The gate
             * above guarantees at least 90% hex characters, so separators are
             * at most sl/10; each separator can cost two entries (a lone
             * nibble and the separator itself), and an entry is a LEB128 index
             * (up to 5 bytes) plus the character — 1.2*sl in the limit.  The
             * old sl/2+16 was not enough: 9 hex digits per line hits the gate
             * exactly and overran it by 3.1MB on an 8MB block. */
            size_t blk = 13 + sl + sl + sl / 4 + 64;
            hexbuf = malloc(blk);
            if (hexbuf) {
                size_t n_sep = 0, bin_len = 0;
                size_t v = 0;                        /* values emitted so far */
                int pend = -1;                       /* pending hex nibble */
                uint32_t orig_len = (uint32_t)sl;
                for (size_t i = 0; i < sl; i++) {
                    uint8_t c = s[i];
                    int hv;
                    if (c >= '0' && c <= '9') hv = c - '0';
                    else if (c >= 'a' && c <= 'f') hv = c - 'a' + 10;
                    else if (c >= 'A' && c <= 'F') hv = c - 'A' + 10;
                    else hv = -1;
                    if (hv >= 0) {
                        if (pend < 0) pend = hv;
                        else { hexbuf[13 + bin_len] = (uint8_t)((pend << 4) | hv); bin_len++; pend = -1; v++; }
                    } else {
                        if (pend >= 0) {             /* lone nibble kept as char */
                            uint32_t p = (uint32_t)v;
                            uint8_t lb[4]; int ln = 0;
                            do { lb[ln++] = (uint8_t)(p & 0x7F); p >>= 7; if (p) lb[ln - 1] |= 0x80; } while (p);
                            for (int b2 = 0; b2 < ln; b2++) hexbuf[13 + sl + n_sep++] = lb[b2];
                            hexbuf[13 + sl + n_sep++] = (uint8_t)(pend < 10 ? '0' + pend : hex_lit + pend - 10);
                            pend = -1;
                        }
                        uint32_t p = (uint32_t)v;
                        uint8_t lb[4]; int ln = 0;
                        do { lb[ln++] = (uint8_t)(p & 0x7F); p >>= 7; if (p) lb[ln - 1] |= 0x80; } while (p);
                        for (int b2 = 0; b2 < ln; b2++) hexbuf[13 + sl + n_sep++] = lb[b2];
                        hexbuf[13 + sl + n_sep++] = c;
                    }
                }
                if (pend >= 0) {
                    uint32_t p = (uint32_t)v;
                    uint8_t lb[4]; int ln = 0;
                    do { lb[ln++] = (uint8_t)(p & 0x7F); p >>= 7; if (p) lb[ln - 1] |= 0x80; } while (p);
                    for (int b2 = 0; b2 < ln; b2++) hexbuf[13 + sl + n_sep++] = lb[b2];
                    hexbuf[13 + sl + n_sep++] = (uint8_t)(pend < 10 ? '0' + pend : hex_lit + pend - 10);
                }
                /* The allocation above is sized for the worst case the gate
                 * can admit; this guards the invariant rather than a known
                 * input, so a future bound change fails loudly instead of
                 * quietly corrupting the heap. */
                if ((size_t)13 + bin_len + n_sep > blk) { free(hexbuf); return -1; }
                memmove(hexbuf + 13 + bin_len, hexbuf + 13 + sl, n_sep);  /* table after bin: lag scan sees clean bin */
                hexbuf[0] = (uint8_t)(orig_len & 0xFF); hexbuf[1] = (uint8_t)((orig_len >> 8) & 0xFF);
                hexbuf[2] = (uint8_t)((orig_len >> 16) & 0xFF); hexbuf[3] = (uint8_t)((orig_len >> 24) & 0xFF);
                hexbuf[4] = hex_case;
                hexbuf[5] = (uint8_t)(bin_len & 0xFF); hexbuf[6] = (uint8_t)((bin_len >> 8) & 0xFF);
                hexbuf[7] = (uint8_t)((bin_len >> 16) & 0xFF); hexbuf[8] = (uint8_t)((bin_len >> 24) & 0xFF);
                hexbuf[9] = (uint8_t)(n_sep & 0xFF); hexbuf[10] = (uint8_t)((n_sep >> 8) & 0xFF);
                hexbuf[11] = (uint8_t)((n_sep >> 16) & 0xFF); hexbuf[12] = (uint8_t)((n_sep >> 24) & 0xFF);
                sl = 13 + bin_len + n_sep;
                s = hexbuf;
                hex_mode = 1;
            }
        }
    }

    /* v5: internal segmentation — mixed blocks (struct tables + text +
     * random...) can't share one lag/stride; a single global pick helps
     * the sampled region and hurts the rest.  Split into 64KB segments,
     * each compressed independently as a sub-block (own lag/stride/tags);
     * segment lengths make the stream self-delimiting.  Segments smaller
     * than the threshold recurse as plain single-block calls.
     *
     * Below 1MB the split used to be skipped, which cost whole blocks: a
     * 384KB byte-plane split of a K tensor (compressible high bytes, then
     * random low bytes) has every group straddling both halves, so no tag
     * can price it and the block fell to RAW at 1.00x where xz reaches
     * 1.17x.  The spread test below segments those without paying the
     * per-segment overhead on blocks that really are homogeneous. */
    int want_seg = (sl >= CPENC_SEG_MIN);
    /* The spread test must not fire on a block that is already one
     * segment: n_seg would be 1 and the recursive call would re-enter
     * this branch with the same size — infinite recursion (a 256KB
     * segment of records overflowed the stack the moment SEG_SIZE was
     * raised above the 128KB test threshold). */
    if (!want_seg && sl > CPENC_SEG_SIZE && sl >= (128u << 10)) {
        size_t w = 1024, nw = 16;
        double lo = 9.0, hi = 0.0;
        if (sl / w >= nw) {
            for (size_t i = 0; i < nw; i++) {
                size_t off = (sl - w) * i / (nw - 1);
                double e = ent(s + off, w);
                if (e < lo) lo = e;
                if (e > hi) hi = e;
            }
            want_seg = (hi - lo) > 1.0;
        }
    }
    /* Segment size.  SEG_SIZE is the content-isolation granularity (big
     * enough that the per-segment table costs fade, small enough to keep
     * mixed regions apart), and it must be a pure function of the DATA:
     * an earlier version widened it to sl/nthr ("one big chunk per
     * core"), which made the encoder's output depend on CPENC_THREADS.
     * The nested calls a top-level block makes — the match layer
     * recompressing its token stream (ref_try), and every segment a
     * worker compresses — run this same code, and ref_try runs on the
     * main thread where the t_in_worker guard is clear, so a thread
     * count of 4 cut an 8MB records block into 2MB chunks at the top and
     * 1MB chunks inside the match layer instead of the 512KB both used
     * at t=1.  4 x 1MB chunks is exactly the pathological split for that
     * token stream: 1MB < 2*SEG_SIZE forces the 64KB fallback, whose
     * chunks then trip the constant-block test (ent < 1.0) and switch
     * segmentation off entirely, and 1373045 bytes became 1540207
     * (6.1095x -> 5.4464x).  The pool splits by the shared counter now,
     * so the per-core chunking buys nothing: 512KB tasks feed 4 workers
     * just as well, and the top level is bit-identical at any thread
     * count. */
    size_t seg_sz = CPENC_SEG_SIZE;
    int nthr = seg_threads(16);
    if (sl < 2 * seg_sz) seg_sz = (64u << 10);   /* too small to split coarsely */
    if (!no_lz) {   /* tuning hook: top level only (nested calls pass no_lz) */
        const char *cs = getenv("CPENC_CHUNK");
        if (cs) { long v = atol(cs); if (v >= 4096) seg_sz = (size_t)v; }
    }
    size_t n_seg = want_seg ? (sl + seg_sz - 1) / seg_sz : 0;
    /* n_seg < 2 would recurse into this same branch at the same size;
     * ent < 1.0 means a constant block, which the normal path encodes in
     * a few bytes. */
    if (want_seg && (n_seg < 2 || ent(s, 2048) < 1.0)) want_seg = 0;
    if (want_seg && (size_t)nthr > n_seg) nthr = (int)n_seg;
    if (want_seg) {
        size_t ext = hex_mode ? 1 : 0;      /* real-stride byte for hex wrapper */
        size_t need = CPENC_HEADER_SIZE + ext + 4 + n_seg * 4;
        if (*dl < need) { free(hexbuf); return -1; }
        hdr_head((uint8_t*)d, hex_mode ? CPENC_STRIDE_HEX : CPENC_STRIDE_SEG, (uint32_t)sl);
        uint8_t *op = d + CPENC_HEADER_SIZE;
        if (hex_mode) *op++ = CPENC_STRIDE_SEG;
        *op++ = (uint8_t)(n_seg & 0xFF); *op++ = (uint8_t)((n_seg >> 8) & 0xFF);
        *op++ = (uint8_t)((n_seg >> 16) & 0xFF); *op++ = (uint8_t)((n_seg >> 24) & 0xFF);
        uint8_t *lenp = op;
        op += n_seg * 4;
        if (nthr > 1) {
            /* parallel: each segment compresses into its own buffer,
             * then a sequential pass assembles the stream */
            uint8_t **bufs = malloc(n_seg * sizeof(uint8_t*));
            SegEncJob *jobs = calloc(n_seg, sizeof(SegEncJob));
            if (!bufs || !jobs) { free(bufs); free(jobs); free(hexbuf); return -1; }
            for (size_t k = 0; k < n_seg; k++) {
                bufs[k] = malloc(cpenc_max_compressed_size(seg_sz) + 64);
                if (!bufs[k]) { for (size_t j = 0; j < k; j++) free(bufs[j]); free(bufs); free(jobs); free(hexbuf); return -1; }
                jobs[k].src = s + k * seg_sz;
                jobs[k].len = sl - k * seg_sz;
                if (jobs[k].len > seg_sz) jobs[k].len = seg_sz;
                jobs[k].out = bufs[k];
                jobs[k].cap = cpenc_max_compressed_size(seg_sz) + 64;
                jobs[k].mo = mo;
            }
            pthread_t th[16];
            SegPool pool;
            pool.src = s; pool.sl = sl; pool.n_seg = n_seg; pool.seg_sz = seg_sz;
            pool.jobs = jobs; pool.next = 0; pool.mo = mo;
            pthread_mutex_init(&pool.mtx, NULL);
            int made = 0;
            for (int i = 0; i < nthr; i++)
                if (pthread_create(&th[i], NULL, seg_pool_worker, &pool) == 0) made++;
            if (!made) {   /* could not start any worker: do it inline */
                pool.next = 0;
                seg_pool_worker(&pool);
            }
            for (int i = 0; i < made; i++) pthread_join(th[i], NULL);
            pthread_mutex_destroy(&pool.mtx);
            for (size_t k = 0; k < n_seg; k++) {
                if (jobs[k].rc != 0 || jobs[k].outlen > (size_t)(*dl - (op - d))) {
                    for (size_t j = 0; j < n_seg; j++) free(bufs[j]);
                    free(bufs); free(jobs); free(hexbuf); return -1;
                }
                lenp[k * 4 + 0] = (uint8_t)(jobs[k].outlen & 0xFF);
                lenp[k * 4 + 1] = (uint8_t)((jobs[k].outlen >> 8) & 0xFF);
                lenp[k * 4 + 2] = (uint8_t)((jobs[k].outlen >> 16) & 0xFF);
                lenp[k * 4 + 3] = (uint8_t)((jobs[k].outlen >> 24) & 0xFF);
                memcpy(op, bufs[k], jobs[k].outlen);
                op += jobs[k].outlen;
            }
            for (size_t j = 0; j < n_seg; j++) free(bufs[j]);
            free(bufs); free(jobs);
        } else {
            for (size_t k = 0; k < n_seg; k++) {
                size_t seg_off = k * seg_sz;
                size_t seg_len = sl - seg_off;
                if (seg_len > seg_sz) seg_len = seg_sz;
                size_t clen = (size_t)(*dl - (op - d));
                if (cpenc_compress_impl(s + seg_off, seg_len, op, &clen, mo, 0, 1) != 0) {
                    free(hexbuf); return -1;
                }
                lenp[k * 4 + 0] = (uint8_t)(clen & 0xFF);
                lenp[k * 4 + 1] = (uint8_t)((clen >> 8) & 0xFF);
                lenp[k * 4 + 2] = (uint8_t)((clen >> 16) & 0xFF);
                lenp[k * 4 + 3] = (uint8_t)((clen >> 24) & 0xFF);
                op += clen;
            }
        }
        hdr_comp((uint8_t*)d, (uint32_t)(op - d - CPENC_HEADER_SIZE));
        *dl = (size_t)(op - d);
        /* Match layer runs on the WHOLE block here, not per segment: a
         * segment has no predecessor, so its first record is literals —
         * for 1KB records that cold start (~1KB) is as large as the rest
         * of the segment's compressed output, halving the ratio (57x -> 31x
         * measured).  Matching across segment boundaries removes it. */
        if (!no_lz && !hex_mode) {
            size_t r = ref_try(s, sl, d, *dl, (size_t)(op - d), mo);
            if (r) *dl = r;
        }
        free(hexbuf);
        return 0;
    }

    /* (v6 BWT candidate moved to the end of this function: it is only
     * worth taking when it actually beats the main path — see the note
     * there.  It must run on the ORIGINAL bytes, so keep them: the lag
     * pass below rebinds `s` to the differenced buffer, and the BWT
     * decode path has no un-diff step.) */
    const uint8_t *s_bwt = s;

    /* v3: lag-delta preprocessing — catches fixed-frame periodicity
     * (log rows, protocol frames, config tables).  Coarse scan (every 16)
     * then fine scan around the best; high-entropy segments skip entirely
     * (no frame periodicity to find).  Diff + entropy merged into one pass. */
    int lag_idx = CPENC_LAG_NONE;
    uint8_t *lagbuf = NULL;
    if (!no_lag && sl >= 4096) {
        size_t k = sl < 8192 ? sl : 8192;
        double best_e = ent(s, k);
        int best_p = 0;
        /* zero-dense streams (BWT+MTF output, sparse data): a lag diff only
         * pollutes — 0 - nonzero wraps to 255 spikes (measured).  Skip the
         * scan entirely; the 0-run machinery handles these. */
        if (best_e > 1.0 && best_e < 6.0) {
            uint32_t zc = 0;
            for (size_t i = 0; i < k; i++) if (!s[i]) zc++;
            if (zc * 2 <= k) {
            /* near-uniform (>6.0) or near-constant (<1.0) segments have no
             * frame periodicity worth scanning for */
            /* coarse: 16 candidates instead of 240, small sample (frame
             * periods 16..255 are visible in 2048 bytes) */
            size_t kc = k < 2048 ? k : 2048;
            for (int p = CPENC_LAG_MIN; p <= CPENC_LAG_MAX && p < (int)kc; p += 16) {
                uint32_t c[256] = {0};
                for (size_t i = 0; i < kc; i++) {
                    uint8_t v = (uint8_t)(s[i] - (i >= (size_t)p ? s[i - (size_t)p] : 0));
                    c[v]++;
                }
                double h = fq_ent(c, kc);
                if (h < best_e - 0.05) { best_e = h; best_p = p; }
            }
            if (best_p) {
                lag_idx = best_p;               /* coarse candidate as floor:
                                                 * fine scan re-tests it at the
                                                 * same entropy, which can't
                                                 * beat the -0.05 margin */
                size_t kf = k < 4096 ? k : 4096;
                int lo = best_p - 8; if (lo < CPENC_LAG_MIN) lo = CPENC_LAG_MIN;
                int hi = best_p + 8; if (hi > CPENC_LAG_MAX) hi = CPENC_LAG_MAX;
                for (int p = lo; p <= hi && p < (int)kf; p++) {
                    uint32_t c[256] = {0};
                    for (size_t i = 0; i < kf; i++) {
                        uint8_t v = (uint8_t)(s[i] - (i >= (size_t)p ? s[i - (size_t)p] : 0));
                        c[v]++;
                    }
                    double h = fq_ent(c, kf);
                    if (h < best_e - 0.05) { best_e = h; lag_idx = p; }
                }
            }
            }
        }
    }
    if (lag_idx) {
        lagbuf = malloc(sl);
        if (!lagbuf) return -2;
        for (size_t i = 0; i < sl; i++)
            lagbuf[i] = (uint8_t)(s[i] - (i >= (size_t)lag_idx ? s[i - (size_t)lag_idx] : 0));
        s = lagbuf;
    }

    /* near-random segments: entropy > 7.5 means no delta/transform can
     * help — store raw.  lag_idx is 0 here (pre-filter skips the scan),
     * so the RAW decode path's lag restore is a no-op. */
    if (sl >= 1024) {
        double seg_e = ent(s, sl);
        /* 7.9, not 7.5: real Q4_K_M weights sit at 7.78 and still have
         * 2.7% of entropy coding in them.  Truly random data is 7.9998, so
         * this still short-circuits the genuinely incompressible case. */
        if (seg_e > 7.9) {
            size_t ext = hex_mode ? 2 : 1;
            if (*dl < CPENC_HEADER_SIZE + ext + sl) { free(lagbuf); free(hexbuf); return -1; }
            /* Near-random *bytes* can still be a repeated record with a few
             * fields changed (byte entropy ~8.0, structure highly
             * redundant) — let the match layer bid before falling to RAW. */
            if(!no_lz && !hex_mode && allow_ref){
                size_t r = ref_try(s_bwt, sl, d, *dl, CPENC_HEADER_SIZE+ext+sl, mo);
                if(r){ *dl=r; free(lagbuf);free(hexbuf); return 0; }
            }
            hdr_put((uint8_t*)d, hex_mode ? CPENC_STRIDE_HEX : CPENC_STRIDE_RAW,
                    (uint32_t)sl,(uint32_t)(ext+sl));
            if(hex_mode) d[CPENC_HEADER_SIZE]=CPENC_STRIDE_RAW;
            d[CPENC_HEADER_SIZE+ext-1]=(uint8_t)lag_idx;
            memcpy(d+CPENC_HEADER_SIZE+ext,s,sl);*dl=CPENC_HEADER_SIZE+ext+sl;
            free(lagbuf);free(hexbuf);return 0;
        }
    }

    int stride=ds(s,sl);
    if(stride<1||stride>8) stride=1;

    /* Stride>1 override: only force stride=1 when raw entropy 3.0-6.0.
     * Above 6.0, the high-byte group in a multi-byte format (int16/float16/
     * float32) often shows strong pattern that stride>1 can exploit. */
    if(stride>1&&sl>=1024){
        double e=ent(s,sl);
        if(e>=3.0&&e<6.0)stride=1;
    }

    if(stride>1&&sl>=256){
        size_t k2 = sl < 1024 ? sl : 1024;
        size_t gs_test = k2/(size_t)stride;
        if(gs_test >= 64){
            uint8_t *tt=malloc(k2),*tb=malloc(gs_test);
            int nraw=0;
            if(tt&&tb){
                tr(s,tt,k2,stride);
                for(int g=0;g<stride;g++){uint8_t *gp=tt+(size_t)g*gs_test;
                    memcpy(tb,gp,gs_test);double e0=ent(tb,gs_test);
                    de(tb,gs_test,1);double e1=ent(tb,gs_test);
                    if(e0>6.5&&e1>6.5)nraw++;}
                if(nraw==stride)stride=1;
            }
            free(tt);free(tb);
        }
    }

    size_t gs=sl/(size_t)stride;
    uint8_t *trans = malloc(sl); if(!trans){free(lagbuf);free(hexbuf);return -2;}
    uint8_t *trial = malloc(4096); if(!trial){free(trans);free(lagbuf);return-2;}

    if(stride==1){memcpy(trans,s,sl);}
    else{tr(s,trans,sl,stride);}

    uint8_t orders[8];
    int nraw=0;

    for(int g=0;g<stride;g++){
        uint8_t *gp = trans + (size_t)g*gs;
        int best_o=0;
        /* One pass builds all four histograms.  Trying each order by
         * copying the sample, differencing it and running an entropy pass
         * costs ~11 traversals of the sample; this loop was 71% of
         * compression time on structured tables.  The chained difference
         * here is the same operator `de` applies in place. */
        size_t sample = gs<4096?gs:4096;
        uint32_t hq[4][256];
        memset(hq, 0, sizeof(hq));
        {
            uint8_t a0p=0, a1p=0, a2p=0;
            for (size_t i2 = 0; i2 < sample; i2++) {
                uint8_t a0 = gp[i2];
                uint8_t a1 = (uint8_t)(a0 - (i2 ? a0p : 0));
                hq[0][a0]++; hq[1][a1]++;
                uint8_t a2 = (uint8_t)(a1 - (i2 ? a1p : 0));
                uint8_t a3 = (uint8_t)(a2 - (i2 ? a2p : 0));
                if (mo>=2) hq[2][a2]++;
                if (mo>=3) hq[3][a3]++;
                a0p=a0; a1p=a1; a2p=a2;
            }
        }
        double sample_e = fq_ent(hq[0], sample);
        double e1 = fq_ent(hq[1], sample);
        if(e1<sample_e-0.3){
            best_o=1; double best_e=e1;
            if(mo>=2 && best_e>1.0){double e2=fq_ent(hq[2],sample);if(e2<best_e-0.01){best_o=2;best_e=e2;}}
            if(mo>=3 && best_e>1.0){double e3=fq_ent(hq[3],sample);if(e3<best_e-0.01){best_o=3;}}
        }else if(sample_e>7.5){best_o=0;nraw++;}   /* 7.5, not 6.5: a 6.6-bit
            * group still has ~17% of entropy coding in it (real FP8 weight
            * blocks sit at 6.67 and were being stored RAW) */
        de(gp,gs,best_o);
        orders[g]=(uint8_t)best_o;
    }

    /* Global RAW only for stride>1 where all groups independently
     * failed delta — a single-group stride=1 with high entropy may
     * still benefit from Huffman/LUT. Let per-group tag selection decide. */
    if(nraw==stride && stride>1 && sl>=256){
        size_t ext = hex_mode ? 2 : 1;
        if(*dl<CPENC_HEADER_SIZE+ext+sl){free(trial);free(trans);free(lagbuf);free(hexbuf);return-1;}
        /* The match layer gets first refusal: RAW is what delta/BWT produce
         * on data that is really a repeated record with a few fields
         * changed, and returning here would hide the match layer's win. */
        if(!no_lz && !hex_mode){
            size_t r = ref_try(s_bwt, sl, d, *dl, CPENC_HEADER_SIZE+ext+sl, mo);
            if(r){ *dl=r; free(trial);free(trans);free(lagbuf);free(hexbuf); return 0; }
        }
        hdr_put((uint8_t*)d, hex_mode ? CPENC_STRIDE_HEX : CPENC_STRIDE_RAW,
                (uint32_t)sl,(uint32_t)(ext+sl));
        if(hex_mode) d[CPENC_HEADER_SIZE]=(uint8_t)stride;
        d[CPENC_HEADER_SIZE+ext-1]=(uint8_t)lag_idx;
        memcpy(d+CPENC_HEADER_SIZE+ext,s,sl);*dl=CPENC_HEADER_SIZE+ext+sl;
        free(trial);free(trans);free(lagbuf);free(hexbuf);return 0;
    }

    uint8_t *out=d+CPENC_HEADER_SIZE,*edd=d+*dl;
    if(out>=edd){free(trial);free(trans);free(lagbuf);free(hexbuf);return-1;}
    if(hex_mode){ if(out>=edd){free(trial);free(trans);free(lagbuf);free(hexbuf);return-1;} *out++=(uint8_t)stride; }   /* real stride for the inner block */
    *out++=(uint8_t)lag_idx;
    for(int g=0;g<stride;g++){if(out>=edd){free(trial);free(trans);free(lagbuf);free(hexbuf);return-1;}*out++=orders[g];}

    for(int g=0;g<stride;g++){
        uint8_t *gp = trans + (size_t)g*gs;
        uint32_t fq[256];
        freq256(gp, gs, fq);
        int nz=0;for(int i=0;i<256;i++)if(fq[i])nz++;

        size_t sz_raw=4+gs; uint8_t tag;
        uint8_t mk_bp[8] = {0};   /* tag 6 plane plan, reused by the emit */
        size_t bcnt_bp[8] = {0};  /* tag 6 per-plane set-bit counts */
        /* Every candidate's bid is its exact emitted byte count, so the
         * tag competition is a true minimum.  Entropy estimates with a
         * safety margin used to forfeit gains under 3% (at 8MB that margin
         * is ~245KB); on real Q4_K_M weights that gave 0.9997x where gzip
         * reached 1.0223x. */
        size_t sz_h = huff_size_exact(gp, gs, fq);

        if(nz<=1){tag=0;}
        else{
            int bits=1;while((1<<bits)<nz&&bits<HUFF_MAX)bits++;
            size_t sz_l=6+(size_t)nz+((size_t)gs*(size_t)bits+7)/8;
            size_t sz_best=sz_raw; tag=3;
            if(sz_l<sz_best){tag=1;sz_best=sz_l;}
            if(sz_h<sz_best){tag=2;sz_best=sz_h;}
            /* tag 6: bit planes.  Packed planes alone can never beat a byte
             * Huffman code (H(byte) <= sum of plane entropies), so the tag
             * only pays where a plane is constant or near-constant.  Plan
             * only when a cheap lower bound can still beat the best bid:
             * per plane, the minority needs at least ceil(mn/8) packed
             * bytes at 2 bytes each, or pn if that is less.  Without this
             * the pack ran on every text group (ASCII has a constant top
             * plane) and halved throughput. */
            {
                size_t pn_bp = (gs + 7) / 8;
                size_t lb = 8;
                for (int b = 0; b < 8; b++) {
                    size_t o = 0;
                    for (int i2 = 0; i2 < 256; i2++)
                        if (fq[i2] && ((i2 >> b) & 1)) o += fq[i2];
                    bcnt_bp[b] = o;
                    if (o == 0 || o == gs) { lb += 1; continue; }
                    size_t mn = o < gs - o ? o : gs - o;
                    size_t sp_lb = 4 + 2 * ((mn + 7) / 8);
                    lb += sp_lb < pn_bp ? sp_lb : pn_bp;
                }
                if (lb < sz_best) {
                    size_t sz_bp = bitplane_plan(gp, gs, bcnt_bp, mk_bp);
                    if (sz_bp != (size_t)-1 && sz_bp < sz_best) {
                        tag = 6; sz_best = sz_bp;
                    }
                }
            }
            /* tag 5: single-value dominance — one symbol covers >99% of the
             * group.  Encode the constant plus a short exception list
             * (positions LEB128 + values), vs Huffman's wasteful 1-bit-per-
             * symbol all-zero bitstream.  Complements RLE (zeros) with
             * non-zero constants (e.g. transposed struct groups). */
            {
                uint32_t dom = 0, domv = 0;
                for (int i2 = 0; i2 < 256; i2++)
                    if (fq[i2] > dom) { dom = fq[i2]; domv = (uint32_t)i2; }
                if (dom > (uint32_t)gs - gs/10 && dom < (uint32_t)gs && gs >= 1024) {
                    /* Exact: [val u16][nex u32] + per exception (gap LEB)(val).
                     * The old 3-byte-per-exception guess both overbid short
                     * gaps and ignored long ones. */
                    size_t sz_c = 6;
                    uint32_t prev = 0; int first = 1;
                    for (size_t i2 = 0; i2 < gs; i2++)
                        if (gp[i2] != (uint8_t)domv) {
                            uint32_t pos = (uint32_t)i2;
                            uint32_t enc = first ? pos : (pos - prev - 1);
                            first = 0; prev = pos;
                            sz_c += (size_t)leb_len(enc) + 1;
                        }
                    if (sz_c < sz_best) { tag = 5; sz_best = sz_c; }
                }
            }
            /* RLE candidate (tag 4): zero-heavy groups get run-length coded.
             * Strict alternation: zero-run then value; consecutive nonzeros
             * are separated by run=0.  Runs are 8-bit with 255 = continuation.
             * Threshold lowered to 50% so BWT+MTF streams (0 ~57-70%, long
             * runs) qualify; run count is scanned (cheap, only for zero-heavy
             * groups) because the alternating-stream bound 2nv+1 overestimates
             * clustered zeros by ~5x and would steer them to Huffman.
             * The decoder sizes its run buffer from this same threshold
             * (rle_max_run_bytes): the encoded run stream is a byte stream
             * with 255+LEB128 escapes, so it is longer than the run count.
             * Loosening gs/2 here means loosening that bound too. */
            if(fq[0] > (uint32_t)gs/2 && gs >= 64){
                uint32_t *runs; uint8_t *vls, *r8; size_t nv, nr8, sz_rle;
                int st_r, st_v;
                if (rle_build(gp, gs, &runs, &vls, &nv, &r8, &nr8, &sz_rle,
                              &st_r, &st_v) == 0) {
                    if (sz_rle < sz_best) { tag = 4; sz_best = sz_rle; }
                    free(r8); free(runs); free(vls);
                }
            }
        }

        if(out>=edd){free(trial);free(trans);free(lagbuf);free(hexbuf);return-1;} *out++=tag;
        if(tag==4){
            /* RLE: the same builder that priced the group; its subtag
             * choices are the ones the bid was computed from. */
            uint32_t *runs; uint8_t *vls, *r8; size_t nv, nr8, sz_rle;
            int st, st_v;
            if (rle_build(gp, gs, &runs, &vls, &nv, &r8, &nr8, &sz_rle,
                          &st, &st_v) != 0) {
                free(trial);free(trans);free(lagbuf);free(hexbuf);return-1;
            }
            uint32_t rq[256]={0};for(size_t k=0;k<nr8;k++)rq[r8[k]]++;
            int rnz=0;for(int i2=0;i2<256;i2++)if(rq[i2])rnz++;
            if(out>=edd){free(r8);free(runs);free(vls);free(trial);free(trans);free(lagbuf);free(hexbuf);return-1;}
            *out++=st;
            if(st==0){
                uint8_t val=0;for(int i2=0;i2<256;i2++)if(rq[i2]){val=(uint8_t)i2;break;}
                if(out+6>edd){free(r8);free(runs);free(vls);free(trial);free(trans);free(lagbuf);free(hexbuf);return-1;}
                *out++=nr8&0xFF;*out++=(nr8>>8)&0xFF;*out++=(nr8>>16)&0xFF;*out++=nr8>>24;
                *out++=val;*out++=(uint8_t)(rnz-1);
            }else if(st==4){
                if(sub5_encode(&out,edd,r8,nr8,rq)!=0){free(r8);free(runs);free(vls);free(trial);free(trans);free(lagbuf);free(hexbuf);return-1;}
            }else if(st==1){
                int rb=1;while((1<<rb)<rnz&&rb<HUFF_MAX)rb++;
                uint8_t lut[256],rev[256]={0};int li=0;
                for(int i2=0;i2<256&&li<rnz;i2++)if(rq[i2]){lut[li]=(uint8_t)i2;rev[i2]=(uint8_t)li;li++;}
                size_t need=6+(size_t)rnz+nr8*(size_t)rb/8+1;
                if((size_t)(edd-out)<need){free(r8);free(runs);free(vls);free(trial);free(trans);free(lagbuf);free(hexbuf);return-1;}
                *out++=nr8&0xFF;*out++=(nr8>>8)&0xFF;*out++=(nr8>>16)&0xFF;*out++=nr8>>24;
                *out++=(uint8_t)(rnz-1);
                for(int i2=0;i2<rnz;i2++)*out++=lut[i2];
                uint64_t enc_buf=0;int enc_bc=0;
                for(size_t k=0;k<nr8;k++){enc_buf|=(uint64_t)rev[r8[k]]<<enc_bc;enc_bc+=rb;
                    while(enc_bc>=8){*out++=(uint8_t)(enc_buf&0xFF);enc_buf>>=8;enc_bc-=8;}}
                if(enc_bc>0)*out++=(uint8_t)(enc_buf&0xFF);
            }else if(st==2){
                uint8_t lens[256];uint16_t code[256];uint8_t bits[256];
                hlens(rq,lens);hcodes(lens,code,bits);
                uint32_t used=h_enc(r8,nr8,lens,code,bits,out,(size_t)(edd-out));
                if(used==0||(size_t)(edd-out)<used){free(r8);free(runs);free(vls);free(trial);free(trans);free(lagbuf);free(hexbuf);return-1;}
                out+=used;
            }else{
                if((size_t)(edd-out)<4+nr8){free(r8);free(runs);free(vls);free(trial);free(trans);free(lagbuf);free(hexbuf);return-1;}
                *out++=nr8&0xFF;*out++=(nr8>>8)&0xFF;*out++=(nr8>>16)&0xFF;*out++=nr8>>24;
                memcpy(out,r8,nr8);out+=nr8;
            }
            /* encode vals sub-block */
            uint32_t vq[256]={0};for(size_t k=0;k<nv;k++)vq[vls[k]]++;
            int vnz=0;for(int i2=0;i2<256;i2++)if(vq[i2])vnz++;
            st=st_v;
            if(out>=edd){free(r8);free(runs);free(vls);free(trial);free(trans);free(lagbuf);free(hexbuf);return-1;}
            *out++=st;
            if(st==0){
                uint8_t val=0;for(int i2=0;i2<256;i2++)if(vq[i2]){val=(uint8_t)i2;break;}
                if(out+6>edd){free(r8);free(runs);free(vls);free(trial);free(trans);free(lagbuf);free(hexbuf);return-1;}
                *out++=nv&0xFF;*out++=(nv>>8)&0xFF;*out++=(nv>>16)&0xFF;*out++=nv>>24;
                *out++=val;*out++=(uint8_t)(vnz-1);
            }else if(st==4){
                if(sub5_encode(&out,edd,vls,nv,vq)!=0){free(r8);free(runs);free(vls);free(trial);free(trans);free(lagbuf);free(hexbuf);return-1;}
            }else if(st==1){
                int vb=1;while((1<<vb)<vnz&&vb<HUFF_MAX)vb++;
                uint8_t lut[256],rev[256]={0};int li=0;
                for(int i2=0;i2<256&&li<vnz;i2++)if(vq[i2]){lut[li]=(uint8_t)i2;rev[i2]=(uint8_t)li;li++;}
                size_t need=6+(size_t)vnz+nv*(size_t)vb/8+1;
                if((size_t)(edd-out)<need){free(r8);free(runs);free(vls);free(trial);free(trans);free(lagbuf);free(hexbuf);return-1;}
                *out++=nv&0xFF;*out++=(nv>>8)&0xFF;*out++=(nv>>16)&0xFF;*out++=nv>>24;
                *out++=(uint8_t)(vnz-1);
                for(int i2=0;i2<vnz;i2++)*out++=lut[i2];
                uint64_t enc_buf=0;int enc_bc=0;
                for(size_t k=0;k<nv;k++){enc_buf|=(uint64_t)rev[vls[k]]<<enc_bc;enc_bc+=vb;
                    while(enc_bc>=8){*out++=(uint8_t)(enc_buf&0xFF);enc_buf>>=8;enc_bc-=8;}}
                if(enc_bc>0)*out++=(uint8_t)(enc_buf&0xFF);
            }else if(st==2){
                uint8_t lens[256];uint16_t code[256];uint8_t bits[256];
                hlens(vq,lens);hcodes(lens,code,bits);
                uint32_t used=h_enc(vls,nv,lens,code,bits,out,(size_t)(edd-out));
                if(used==0||(size_t)(edd-out)<used){free(r8);free(runs);free(vls);free(trial);free(trans);free(lagbuf);free(hexbuf);return-1;}
                out+=used;
            }else{
                if((size_t)(edd-out)<4+nv){free(r8);free(runs);free(vls);free(trial);free(trans);free(lagbuf);free(hexbuf);return-1;}
                *out++=nv&0xFF;*out++=(nv>>8)&0xFF;*out++=(nv>>16)&0xFF;*out++=nv>>24;
                memcpy(out,vls,nv);out+=nv;
            }
            free(r8);free(runs);free(vls);
        }else if(tag==0){
            uint8_t val=0;for(int i=0;i<256;i++)if(fq[i]){val=(uint8_t)i;break;}
            if(out+6>edd){free(trial);free(trans);free(lagbuf);free(hexbuf);return-1;}
            *out++=gs&0xFF;*out++=(gs>>8)&0xFF;*out++=(gs>>16)&0xFF;*out++=gs>>24;
            *out++=val;*out++=(uint8_t)(nz-1);
        }else if(tag==5){
            /* const + exceptions: [val][nex 4B][pos LEB128 stream][val stream] */
            uint32_t domv2 = 0, dom2 = 0;
            for (int i2 = 0; i2 < 256; i2++) if (fq[i2] > dom2) { dom2 = fq[i2]; domv2 = (uint32_t)i2; }
            uint32_t nex = (uint32_t)gs - dom2;
            if (out+6 > edd){free(trial);free(trans);free(lagbuf);free(hexbuf);return-1;}
            *out++=(uint8_t)(domv2&0xFF);*out++=(uint8_t)(domv2>>8);
            *out++=nex&0xFF;*out++=(nex>>8)&0xFF;*out++=(nex>>16)&0xFF;*out++=nex>>24;
            /* exception positions: first absolute, then gaps (avg gap is
             * small for zero-dominated data — 1 LEB128 byte vs 3-4 for
             * absolute) */
            uint32_t prev = 0;
            int first = 1;
            for (size_t i2 = 0; i2 < gs; i2++) {
                if (gp[i2] != (uint8_t)domv2) {
                    uint32_t pos = (uint32_t)i2;
                    uint32_t enc = first ? pos : (pos - prev - 1);
                    first = 0; prev = pos;
                    uint8_t lb[4]; int ln = 0;
                    do { lb[ln++] = (uint8_t)(enc & 0x7F); enc >>= 7;
                         if (enc) lb[ln-1] |= 0x80; } while (enc);
                    for (int b2 = 0; b2 < ln; b2++) { if (out >= edd){free(trial);free(trans);free(lagbuf);free(hexbuf);return-1;} *out++ = lb[b2]; }
                    if (out >= edd){free(trial);free(trans);free(lagbuf);free(hexbuf);return-1;}
                    *out++ = gp[i2];
                }
            }
        }else if(tag==1){
            int bits=1;while((1<<bits)<nz&&bits<HUFF_MAX)bits++;
            uint8_t lut[256],rev[256]={0};int li=0;
            for(int i=0;i<256&&li<nz;i++)if(fq[i]){uint8_t v=(uint8_t)i;lut[li]=v;rev[v]=(uint8_t)li;li++;}
            size_t need=6+(size_t)nz+gs*(size_t)bits/8+1;
            if((size_t)(edd-out)<need){free(trial);free(trans);free(lagbuf);free(hexbuf);return-1;}
            *out++=gs&0xFF;*out++=(gs>>8)&0xFF;*out++=(gs>>16)&0xFF;*out++=gs>>24;
            *out++=(uint8_t)(nz-1);
            for(int i=0;i<nz;i++)*out++=lut[i];
            /* Wide scalar LUT encode — 8 bytes per batch */
            uint64_t enc_buf=0;int enc_bc=0;
            size_t enc_i;
            for(enc_i=0;enc_i+8<=gs;enc_i+=8){
                uint64_t ei0=rev[gp[enc_i]],ei1=rev[gp[enc_i+1]],ei2=rev[gp[enc_i+2]],ei3=rev[gp[enc_i+3]],
                         ei4=rev[gp[enc_i+4]],ei5=rev[gp[enc_i+5]],ei6=rev[gp[enc_i+6]],ei7=rev[gp[enc_i+7]];
                enc_buf|=ei0<<enc_bc;enc_bc+=bits;while(enc_bc>=8){*out++=(uint8_t)(enc_buf&0xFF);enc_buf>>=8;enc_bc-=8;}
                enc_buf|=ei1<<enc_bc;enc_bc+=bits;while(enc_bc>=8){*out++=(uint8_t)(enc_buf&0xFF);enc_buf>>=8;enc_bc-=8;}
                enc_buf|=ei2<<enc_bc;enc_bc+=bits;while(enc_bc>=8){*out++=(uint8_t)(enc_buf&0xFF);enc_buf>>=8;enc_bc-=8;}
                enc_buf|=ei3<<enc_bc;enc_bc+=bits;while(enc_bc>=8){*out++=(uint8_t)(enc_buf&0xFF);enc_buf>>=8;enc_bc-=8;}
                enc_buf|=ei4<<enc_bc;enc_bc+=bits;while(enc_bc>=8){*out++=(uint8_t)(enc_buf&0xFF);enc_buf>>=8;enc_bc-=8;}
                enc_buf|=ei5<<enc_bc;enc_bc+=bits;while(enc_bc>=8){*out++=(uint8_t)(enc_buf&0xFF);enc_buf>>=8;enc_bc-=8;}
                enc_buf|=ei6<<enc_bc;enc_bc+=bits;while(enc_bc>=8){*out++=(uint8_t)(enc_buf&0xFF);enc_buf>>=8;enc_bc-=8;}
                enc_buf|=ei7<<enc_bc;enc_bc+=bits;while(enc_bc>=8){*out++=(uint8_t)(enc_buf&0xFF);enc_buf>>=8;enc_bc-=8;}
            }
            for(;enc_i<gs;enc_i++){enc_buf|=(uint64_t)rev[gp[enc_i]]<<enc_bc;enc_bc+=bits;
                while(enc_bc>=8){*out++=(uint8_t)(enc_buf&0xFF);enc_buf>>=8;enc_bc-=8;}}
            if(enc_bc>0)*out++=(uint8_t)(enc_buf&0xFF);
        }else if(tag==2){
            uint32_t used;
            if(gs >= (1u<<20))
                used = h_enc_seg(gp, gs, out, (size_t)(edd-out));
            else {
                uint8_t lens[256];uint16_t code[256];uint8_t bits[256];
                hlens(fq,lens);hcodes(lens,code,bits);
                used = h_enc(gp,gs,lens,code,bits,out,(size_t)(edd-out));
            }
            if(used==0||(size_t)(edd-out)<used){free(trial);free(trans);free(lagbuf);free(hexbuf);return-1;}
            out+=used;
        }else if(tag==6){
            /* bit planes — markers were planned during selection */
            size_t pn_bp = (gs + 7) / 8;
            if ((size_t)(edd - out) < 16) { free(trial);free(trans);free(lagbuf);free(hexbuf); return -1; }
            uint8_t *mk = out; out += 8;
            memcpy(mk, mk_bp, 8);
            uint8_t *pk = malloc(pn_bp ? pn_bp : 1);
            if (!pk) { free(trial);free(trans);free(lagbuf);free(hexbuf); return -1; }
            for (int b = 0; b < 8; b++) {
                if (mk[b] == 0) { *out++ = (uint8_t)(bcnt_bp[b] ? 1 : 0); continue; }
                plane_pack(gp, gs, b, pk);
                if (mk[b] == 3)
                    for (size_t i2 = 0; i2 < pn_bp; i2++) pk[i2] ^= 0xFF;
                if (mk[b] >= 2) {
                    size_t nz = 0;
                    for (size_t i2 = 0; i2 < pn_bp; i2++) if (pk[i2]) nz++;
                    if ((size_t)(edd - out) < 4 + nz*6) {
                        free(pk);free(trial);free(trans);free(lagbuf);free(hexbuf); return -1;
                    }
                    *out++ = (uint8_t)(nz & 0xFF); *out++ = (uint8_t)((nz >> 8) & 0xFF);
                    *out++ = (uint8_t)((nz >> 16) & 0xFF); *out++ = (uint8_t)((nz >> 24) & 0xFF);
                    uint32_t prev = 0;
                    for (size_t i2 = 0; i2 < pn_bp; i2++) if (pk[i2]) {
                        out += leb_put(out, (uint32_t)i2 - prev);
                        *out++ = pk[i2]; prev = (uint32_t)i2;
                    }
                } else {
                    if ((size_t)(edd - out) < pn_bp) {
                        free(pk);free(trial);free(trans);free(lagbuf);free(hexbuf); return -1;
                    }
                    memcpy(out, pk, pn_bp); out += pn_bp;
                }
            }
            free(pk);
        }else{
            if((size_t)(edd-out)<4+gs){free(trial);free(trans);free(lagbuf);free(hexbuf);return-1;}
            *out++=gs&0xFF;*out++=(gs>>8)&0xFF;*out++=(gs>>16)&0xFF;*out++=gs>>24;
            memcpy(out,gp,gs);out+=gs;
        }
    }

    hdr_head((uint8_t*)d, hex_mode ? CPENC_STRIDE_HEX : (uint8_t)stride, (uint32_t)sl);
    /* v2: store the stride remainder (tail bytes not covered by any group)
     * so decompress can restore them exactly */
    {
        size_t gend = gs * (size_t)stride;
        if (gend < sl) {
            size_t rem = sl - gend;
            if ((size_t)(edd - out) < rem) { free(trial); free(trans); free(hexbuf); return -1; }
            memcpy(out, trans + gend, rem);
            out += rem;
        }
    }
    /* v7: match-layer candidate — exact repeats as back-references.  Logs
     * are a fixed template repeated per record, which the delta/RLE family
     * cannot express: measured xz 116x vs cpenc 2.9x on a repetitive log,
     * while cpenc wins 317x vs 24x on structured tables.  Built after the
     * main path and kept only if smaller, so it can never make a block
     * worse; the main_len gate skips it where it loses (structured tables,
     * source code and access logs all beat it).
     * Format: [payload], payload = compressed token stream. */
    int lz_won = 0;
    if (!no_lz && !hex_mode && allow_ref) {
        /* Do NOT return on a win: the BWT candidate below still gets to
         * compete against this result.  On source code the reference layer lands between
         * the main path and BWT, and returning here hid BWT's win
         * (4.65x -> 3.30x). */
        size_t r = ref_try(s_bwt, sl, d, *dl, (size_t)(out - d), mo);
        if (r) { out = d + r; lz_won = 1; }
    }

    /* v6: BWT candidate — reorder so context-similar bytes are adjacent,
     * then the lag/stride machinery eats the order-1 structure ("BWT is
     * the delta aligner").  Built AFTER the main path and kept only if it
     * is strictly smaller: a byte-entropy bound is a poor proxy for what
     * delta+transpose already achieves (a structured u32 table has raw
     * byte entropy ~7.1 but collapses to a few hundred bytes), and
     * accepting on that bound lost 317x -> 68x.  Archive-oriented: the
     * sort is slow, so it must earn its place.
     *
     * v7 (stride 0xF9) adds a second candidate.  Format:
     *   [u8 k][eof_row u32][k * anchor u32][inner cpenc block]
     * where the inner block is the RLE2-coded MTF stream, compressed with
     * no_bwt=1, no_lz=1 so it never re-enters this candidate.  The anchors
     * let decode run k LF chains at once.  It competes with the v6 0xFC
     * stream (MTF only, serial decode) and the smaller block is emitted;
     * 0xFC therefore still appears, and still decodes through its original
     * path. */
    /* Probe measured unconditionally under CPENC_VFEAT so the REF-winning
     * blocks get a drop value too — this is the experiment that decides
     * whether the drop can act as a pre-REF discriminator. */
    if (getenv("CPENC_VFEAT") && !no_bwt && !hex_mode && sl >= 8192 && sl <= (1u<<20)) {
        size_t pn = sl < 4096 ? sl : 4096;
        uint8_t *pb = malloc(pn + 1);
        if (pb) {
            bwt_enc(s, pn, pb);
            mtf_enc_inplace(pb, pn + 1);
            uint32_t pc[256] = {0};
            for (size_t i = 0; i <= pn; i++) pc[pb[i]]++;
            fprintf(stderr, "PROBE idx=%s n=%zu drop=%.4f\n", getenv("CPENC_BLKIDX") ? getenv("CPENC_BLKIDX") : "-", sl, ent(s, pn) - fq_ent(pc, pn + 1));
            free(pb);
        }
    }
    if (!no_bwt && !hex_mode && allow_bwt && sl >= 8192 && sl <= (1u<<20)) {
        double e = ent(s_bwt, 2048);
        size_t main_len = (size_t)(out - d);
        /* Two cheap gates before paying for the sort (~4ms per 64KB):
         *  - main_len*8 > sl: BWT tops out around 4-5x, so when the main
         *    path already beats 8x the sort cannot pay for itself (keeps
         *    structured tables at 317x instead of sorting to reject).
         *  - probe: sort a 4KB sample and require MTF to actually
         *    collapse it.  f32-style data passes the entropy window on
         *    ~248 of 256 segments and loses every one of them; the probe
         *    rejects those at 1/16th the sort cost (measured MTF drop:
         *    text +2.85, f32 -0.12 — a clean separation the entropy
         *    window cannot make, since their per-segment byte entropies
         *    overlap at 5.1-6.1).
         * The old lower bound e > 4.7 is gone: a block whose first 2KB are
         * very repetitive but whose body is not (a header over a log, mixed
         * segments) was being skipped, and with anchors+RLE2 the BWT
         * candidate now beats the main path there.  Measured on the 16MB
         * text corpus: 3 of 32 blocks flip to BWT2, 5.52x -> 5.60x, while
         * floats/random/binary are unchanged (their prefix entropy is above
         * the old bound anyway and the probe still rejects them).  The
         * upper bound stays: above 7.2 bits/byte there is nothing to sort. */
        if (e < 7.2 && main_len > CPENC_HEADER_SIZE + 4 + 64 &&
            main_len * 8 > sl) {
            size_t pn = sl < 4096 ? sl : 4096;
            uint8_t *pb = malloc(pn + 1);
            int probe_ok = 0;
            if (pb) {
                bwt_enc(s_bwt, pn, pb);
                mtf_enc_inplace(pb, pn + 1);
                uint32_t pc[256] = {0};
                for (size_t i = 0; i <= pn; i++) pc[pb[i]]++;
                double e0 = ent(s_bwt, pn);
                double e1 = fq_ent(pc, pn + 1);
                probe_ok = (e1 < e0 - 0.3);
                if (getenv("CPENC_VFEAT"))
                    fprintf(stderr, "PROBE n=%zu drop=%.4f ok=%d\n", sl, e0 - e1, probe_ok);
                free(pb);
            }
            uint8_t *bw = probe_ok ? malloc(sl + 1) : NULL;
            if (bw) {
                size_t oi = bwt_enc(s_bwt, sl, bw);
                /* v7: record k anchors along the LF cycle *before* MTF
                 * destroys the L column.  anchor[j] is the position j*m
                 * steps along the cycle (m = (sl+1)/k); the decoder starts
                 * one chain per anchor, which breaks the serial latency of
                 * the single LF chain.  8 anchors already reach the
                 * memory-level-parallelism plateau; 16 help on 1MB blocks. */
                int k = (sl <= (512u << 10)) ? 8 : 16;
                uint32_t anc[16];
                int anc_ok = bwt_make_anchors(bw, sl + 1, oi, anc, k) == 0;
                /* bzip2-style: MTF on the BWT stream turns its local
                 * clustering into long zero runs (RLE territory) */
                mtf_enc_inplace(bw, sl + 1);
                /* RLE2: fold those zero runs before the Huffman pass.
                 * MTF ranks 254/255 cost two bytes each, everything else
                 * stays one, so this is the bzip2 scheme. */
                uint8_t *rle = anc_ok ? malloc(2 * (sl + 1) + 8) : NULL;
                size_t rn = rle ? rle2_enc(bw, sl + 1, rle, 2 * (sl + 1) + 8) : 0;
                /* Never trust a transform blind: decode the run stream back
                 * and compare before spending a Huffman pass on it.  The
                 * standalone prototype's self-test caught a real bug, so it
                 * ships in the encoder too; a mismatch drops the candidate
                 * rather than emitting a stream that cannot round-trip. */
                int rle_ok = 0;
                if (rn) {
                    uint8_t *chk = malloc(sl + 1);
                    if (chk) {
                        rle_ok = (rle2_dec(rle, rn, chk, sl + 1) == sl + 1 &&
                                  memcmp(chk, bw, sl + 1) == 0);
                        free(chk);
                    }
                }
                if (rle_ok) {
                    /* One candidate, chosen by the only cheap signal that
                     * predicts the winner: whether RLE2 actually shortens the
                     * MTF stream.  MTF ranks 254/255 cost two bytes each, so
                     * RLE2 can expand; when it does, the plain MTF stream
                     * wins.  Measured on 16MB of text: RLE2 wins 30 of 32
                     * blocks, records 16 of 16, and the two MTF wins are
                     * 0.15% and 0.66% of their block.  Compressing both
                     * candidates bought 0.02% ratio for 2x encode, so the
                     * loser is not computed any more. */
                    int use_rle = (rn < sl + 1);
                    size_t ovh = use_rle ? (1 + 4 + (size_t)k * 4) : 4;
                    size_t bcap = cpenc_max_compressed_size(use_rle ? rn : sl + 1) + 64;
                    uint8_t *bb = malloc(bcap);
                    if (bb) {
                        size_t cl = bcap - (CPENC_HEADER_SIZE + ovh);
                        int rc;
                        if (use_rle)
                            rc = cpenc_compress_core(rle, rn, bb + CPENC_HEADER_SIZE + ovh,
                                                     &cl, mo, 1, 1, 1);
                        else
                            rc = cpenc_compress_impl(bw, sl + 1, bb + CPENC_HEADER_SIZE + ovh,
                                                     &cl, mo, 1, 1);
                        if (rc == 0 && CPENC_HEADER_SIZE + ovh + cl < main_len) {
                            uint8_t *p = bb + CPENC_HEADER_SIZE;
                            if (use_rle) {
                                hdr_put((uint8_t*)bb, CPENC_STRIDE_BWT2, (uint32_t)sl,
                                        (uint32_t)(ovh + cl));
                                *p++ = (uint8_t)k;
                                *p++ = (uint8_t)oi;
                                *p++ = (uint8_t)(oi >> 8);
                                *p++ = (uint8_t)(oi >> 16);
                                *p++ = (uint8_t)(oi >> 24);
                                for (int j = 0; j < k; j++) {
                                    *p++ = (uint8_t)anc[j];
                                    *p++ = (uint8_t)(anc[j] >> 8);
                                    *p++ = (uint8_t)(anc[j] >> 16);
                                    *p++ = (uint8_t)(anc[j] >> 24);
                                }
                            } else {
                                hdr_put((uint8_t*)bb, CPENC_STRIDE_BWT, (uint32_t)sl,
                                        (uint32_t)(ovh + cl));
                                *p++ = (uint8_t)oi;
                                *p++ = (uint8_t)(oi >> 8);
                                *p++ = (uint8_t)(oi >> 16);
                                *p++ = (uint8_t)(oi >> 24);
                            }
                            memcpy(d, bb, CPENC_HEADER_SIZE + ovh + cl);
                            *dl = CPENC_HEADER_SIZE + ovh + cl;
                            free(bb); free(rle); free(bw);
                            free(trial); free(trans); free(lagbuf); free(hexbuf);
                            return 0;
                        }
                        free(bb);
                    }
                }
                free(rle);
                free(bw);
            }
        }
    }

    /* The reference candidate wrote its own header into d; BWT (which returns
     * directly on a win) had first refusal on beating it.  Nothing left to
     * do but publish the size — rewriting the header here would clobber
     * the candidate's stride. */
    if (lz_won) {
        *dl = (size_t)(out - d);
        free(trial); free(trans); free(lagbuf); free(hexbuf);
        return 0;
    }

    hdr_head((uint8_t*)d, hex_mode ? CPENC_STRIDE_HEX : (uint8_t)stride, (uint32_t)sl);
    hdr_comp((uint8_t*)d, (uint32_t)(out-(d+CPENC_HEADER_SIZE)));
    *dl=(size_t)(out-d); free(trial); free(trans); free(lagbuf); free(hexbuf); return 0;
}

/* Upper bound on the ENCODED run stream of a tag-4 (RLE) group — which is
 * not the run count.  rle_build() writes a run of 255 or more as the byte 255
 * followed by a LEB128 remainder, so one run can occupy several bytes; the
 * decoder walks those bytes, so what it has to size its buffer for is the
 * byte count, not the number of runs.
 *
 * The run count itself is bounded by the encoder's zero-density gate
 * (fq[0] > gs/2, in the RLE candidate of cpenc_compress_impl).  Runs and
 * values strictly alternate, so with z zero-blocks, b non-zero blocks and nv
 * non-zero bytes the run count is z + (nv-b) + [gp[0]!=0]; a group starting
 * with a zero block has z <= b+1 and one starting with a non-zero has z <= b,
 * so either way nr <= nv+1.  The gate caps nv at ceil(gs/2)-1.
 * Each run of 255 or more then adds at most 5 continuation bytes and eats at
 * least 255 zeros, so at most gs/255 of them exist.
 *
 * Both halves come from that one gate: if it is ever loosened, this bound
 * has to be loosened with it.  (It was loosened once — to gs/8 — and the
 * encoder's own runs[] buffer is what overflowed; see the note in
 * rle_build.)  The bound is a corrupt-input guard, so erring large is safe. */
static size_t rle_max_run_bytes(size_t gs) {
    return gs/2 + 2 + 5*(gs/255) + 8;
}

/* ── Nested-call range check ─────────────────────────────────────
 * Every nested decompress works on bytes its caller was given: either a
 * slice of that stream, or a fresh copy (the hex layer makes one).  A call
 * that is neither has had its pointer or its length computed from an
 * untrusted field incorrectly — and once that happens every bounds check
 * further down is arithmetic on a lie, which is how a corrupted block
 * reached h_dec with an end pointer beyond its buffer.  Checking here costs
 * a few compares per block, not per byte, and refuses the block instead. */
#define CPENC_MAX_NEST 8
static __thread const void *g_rb[CPENC_MAX_NEST];
static __thread size_t      g_rn[CPENC_MAX_NEST];
static __thread int         g_rd;

static int cpenc_decompress_impl(const void *sv, size_t sl, void *dv, size_t *dl);

int cpenc_decompress(const void *sv, size_t sl, void *dv, size_t *dl) {
    uintptr_t b0 = (uintptr_t)sv, b1 = b0 + sl;
    /* Nesting has a natural bound: hex and the reference layers each wrap
     * one inner block, so a real block is at most two or three deep.  A
     * crafted one can nest arbitrarily far and blow the stack (fuzzing found
     * that as a stack-overflow before it found anything else). */
    if (g_rd >= CPENC_MAX_NEST) return -1;
    if (g_rd > 0) {
        uintptr_t p0 = (uintptr_t)g_rb[g_rd - 1], p1 = p0 + g_rn[g_rd - 1];
        int inside   = (b0 >= p0 && b1 <= p1);
        int disjoint = (b1 <= p0 || b0 >= p1);
        if (!inside && !disjoint) {
#ifdef CPENC_RANGE_REPORT
            fprintf(stderr, "RANGE: [%p,+%zu) neither within nor disjoint from [%p,+%zu) ra=%p\n",
                    (void *)b0, sl, (void *)p0, g_rn[g_rd - 1], __builtin_return_address(0));
            abort();
#endif
            return -1;
        }
    }
    int save = g_rd;
    if (g_rd < CPENC_MAX_NEST) { g_rb[g_rd] = sv; g_rn[g_rd] = sl; }
    g_rd++;
    int rc = cpenc_decompress_impl(sv, sl, dv, dl);
    g_rd = save;
    return rc;
}

static int cpenc_decompress_impl(const void *sv, size_t sl, void *dv, size_t *dl) {
    const uint8_t *s=(const uint8_t*)sv; uint8_t *d=(uint8_t*)dv;
    if(sl<CPENC_HEADER_SIZE) return -1;
    /* Copy the header out rather than casting: the API takes const void*, so
     * a caller may legitimately hand over a slice that is not 4-byte
     * aligned, and the cast made that undefined behaviour (fuzzing found it
     * with UBSan).  The container always passes malloc'd blocks, but there
     * is no reason for the contract to require it. */
    CPEncHeader hb;
    memcpy(&hb, s, sizeof hb);
    const CPEncHeader *h = &hb;
    /* v1 blocks (no stride remainder) remain decodable; reject only
     * unknown versions and garbage magic */
    if(h->magic!=CPENC_MAGIC||h->version<1||h->version>CPENC_VERSION) return -1;
    if(h->original_size > CPENC_MAX_BLOCK) return -1;
    /* Same ceiling on the compressed side.  A block cannot legitimately be
     * larger than the data it holds (RAW costs original_size + 13), so
     * anything past this is a caller handing over a length it computed from
     * a corrupt field — and every subsequent bounds check inside would be
     * arithmetic on a pointer that far out of range. */
    if(sl > (size_t)CPENC_MAX_BLOCK + 64) return -1;
    /* Every stride except HEX and SEG produces exactly h->original_size bytes,
     * so the caller's capacity can be checked once, here, rather than trusting
     * each branch to remember.  BWT2 (0xF9) is in that flat group: it emits
     * exactly original_size bytes.  It did not: the BWT branch wrote
     * original_size bytes into whatever buffer it was handed, which is a heap
     * overflow when the caller sized that buffer from a different (equally
     * untrusted) number.  HEX and SEG are exempt because they legitimately
     * expand beyond the declared size — a hex sub-block restores to about
     * twice its folded length — and carry their own checks. */
    if (h->stride != CPENC_STRIDE_HEX && h->stride != CPENC_STRIDE_SEG &&
        *dl < (size_t)h->original_size)
        return -1;
    if(h->stride==CPENC_STRIDE_SEG){
        /* v5: independent sub-blocks, lengths in a u32 array */
        const uint8_t *p = s + CPENC_HEADER_SIZE;
        if(sl < CPENC_HEADER_SIZE + 4) return -1;
        uint32_t n_seg = (uint32_t)p[0] | ((uint32_t)p[1]<<8) | ((uint32_t)p[2]<<16) | ((uint32_t)p[3]<<24);
        p += 4;
        /* The encoder splits at CPENC_SEG_SIZE (512KB), except below 1MB,
         * where it drops to 64KB to keep mixed blocks apart — and CPENC_CHUNK
         * may lower that to 4KB.  The bound has to be the floor, not the
         * default: bounding it by 512KB rejected real encoder output (a
         * 512KB REF2 column arrives as 9 segments, not 1) and made blocks
         * that had compressed fine impossible to read back.  The old /4
         * bound was the opposite mistake — 64M segments in a 256MB block,
         * multi-gigabyte job arrays (fuzzing found that at 1.9GB RSS). */
        uint32_t max_seg = h->original_size / (4u << 10) + 2;
        if(n_seg < 1 || n_seg > max_seg) return -1;
        if(sl < CPENC_HEADER_SIZE + 4 + (size_t)n_seg * 4) return -1;
        const uint8_t *sp = p + (size_t)n_seg * 4;
        /* validate all lengths up front */
        for(uint32_t k = 0; k < n_seg; k++) {
            uint32_t clen = (uint32_t)p[k*4] | ((uint32_t)p[k*4+1]<<8) | ((uint32_t)p[k*4+2]<<16) | ((uint32_t)p[k*4+3]<<24);
            if((size_t)(sp + clen) > (size_t)s + sl) return -1;
            sp += clen;
        }
        sp = p + (size_t)n_seg * 4;
        int nthr = seg_threads(n_seg);
        if (nthr > 1) {
            /* parallel decode into per-segment buffers, then assemble */
            uint8_t **bufs = malloc(n_seg * sizeof(uint8_t*));
            SegDecJob *jobs = calloc(n_seg, sizeof(SegDecJob));
            size_t *wants = calloc(n_seg, sizeof(size_t));
            if (!bufs || !jobs || !wants) { free(bufs); free(jobs); free(wants); return -1; }
            for (uint32_t k = 0; k < n_seg; k++) {
                uint32_t clen = (uint32_t)p[k*4] | ((uint32_t)p[k*4+1]<<8) | ((uint32_t)p[k*4+2]<<16) | ((uint32_t)p[k*4+3]<<24);
                /* Segment sizes are an encoder choice; the stream only
                 * stores each sub-block's compressed length.  Its own
                 * header carries original_size, so read it rather than
                 * assuming a fixed granularity. */
                if (clen < CPENC_HEADER_SIZE) { for (uint32_t j = 0; j < k; j++) free(bufs[j]); free(bufs); free(jobs); free(wants); return -1; }
                size_t dl_k = (size_t)sp[4] | ((size_t)sp[5] << 8) | ((size_t)sp[6] << 16) | ((size_t)sp[7] << 24);
                size_t cap = seg_out_cap(sp, dl_k);
                wants[k] = dl_k;
                bufs[k] = malloc(cap ? cap : 1);
                if (!bufs[k]) { for (uint32_t j = 0; j < k; j++) free(bufs[j]); free(bufs); free(jobs); free(wants); return -1; }
                jobs[k].src = sp; jobs[k].clen = clen;
                jobs[k].out = bufs[k]; jobs[k].cap = cap;
                sp += clen;
            }
            pthread_t th[16];
            for (uint32_t base = 0; base < n_seg; base += (uint32_t)nthr) {
                int nb = (int)(n_seg - base < (uint32_t)nthr ? n_seg - base : (uint32_t)nthr);
                for (int i = 0; i < nb; i++) pthread_create(&th[i], NULL, seg_dec_worker, &jobs[base + i]);
                for (int i = 0; i < nb; i++) pthread_join(th[i], NULL);
            }
            size_t o = 0;
            for (uint32_t k = 0; k < n_seg; k++) {
                if (jobs[k].rc != 0 || jobs[k].outlen > jobs[k].cap ||
                    (jobs[k].src[3] != CPENC_STRIDE_HEX && jobs[k].outlen != wants[k])) {
                    for (uint32_t j = 0; j < n_seg; j++) free(bufs[j]);
                    free(bufs); free(jobs); free(wants); return -1;
                }
                if ((size_t)*dl < o + jobs[k].outlen) {
                    for (uint32_t j = 0; j < n_seg; j++) free(bufs[j]);
                    free(bufs); free(jobs); free(wants); return -1;
                }
                memcpy(d + o, bufs[k], jobs[k].outlen);
                o += jobs[k].outlen;
            }
            for (uint32_t j = 0; j < n_seg; j++) free(bufs[j]);
            free(bufs); free(jobs); free(wants);
            *dl = o;
            return 0;
        }
        size_t o = 0;
        for(uint32_t k = 0; k < n_seg; k++) {
            uint32_t clen = (uint32_t)p[k*4] | ((uint32_t)p[k*4+1]<<8) | ((uint32_t)p[k*4+2]<<16) | ((uint32_t)p[k*4+3]<<24);
            /* The parallel path rejects a segment shorter than a header before
             * reading it; this one read four bytes of sp[4..7] regardless, off
             * the end of a short segment near the block's tail. */
            if (clen < CPENC_HEADER_SIZE) return -1;
            size_t dl_k = (size_t)sp[4] | ((size_t)sp[5] << 8) | ((size_t)sp[6] << 16) | ((size_t)sp[7] << 24);
            /* decode into the remaining output space: a hex-wrapped segment
             * restores to more than its declared size */
            size_t cap = (size_t)*dl - o;
            if (cap < dl_k) return -1;
            size_t got = cap;
            if(cpenc_decompress(sp, clen, d + o, &got) != 0) return -1;
            if (got < dl_k || got > cap) return -1;
            if (sp[3] != CPENC_STRIDE_HEX && got != dl_k) return -1;
            o += got;
            sp += clen;
        }
        *dl = o;
        return 0;
    }
    if(h->stride==CPENC_STRIDE_HEX){
        /* v4 hex layer: the payload is a standard block compressing
         * [bin_len u32][n_sep u32][(idx LEB128)(char)]*n_sep [bin bytes].
         * Decompress it as an ordinary block (stride forced to 1 so this
         * branch isn't re-entered), then restore the hex-ascii layout. */
        size_t blk = h->original_size;
        if(sl < CPENC_HEADER_SIZE + 1) return -1;
        int rs = s[CPENC_HEADER_SIZE];      /* real stride of the inner block */
        if(!((rs >= 1 && rs <= 8) || rs == CPENC_STRIDE_SEG || rs == CPENC_STRIDE_RAW)) return -1;
        uint8_t *ns = malloc(sl);
        if(!ns) return -1;
        memcpy(ns, s, sl);
        ((CPEncHeader*)ns)->stride = (uint8_t)rs;
        memmove(ns + CPENC_HEADER_SIZE, ns + CPENC_HEADER_SIZE + 1, sl - CPENC_HEADER_SIZE - 1);
        uint8_t *tmp = malloc(blk ? blk : 1);
        if(!tmp){ free(ns); return -1; }
        size_t tdl = blk;
        if(cpenc_decompress(ns, sl, tmp, &tdl) != 0){ free(ns); free(tmp); return -1; }
        free(ns);
        if(tdl < 13){ free(tmp); return -1; }
        uint32_t orig_len = (uint32_t)tmp[0] | ((uint32_t)tmp[1]<<8) | ((uint32_t)tmp[2]<<16) | ((uint32_t)tmp[3]<<24);
        int hex_case = tmp[4] ? 1 : 0;
        uint32_t bin_len = (uint32_t)tmp[5] | ((uint32_t)tmp[6]<<8) | ((uint32_t)tmp[7]<<16) | ((uint32_t)tmp[8]<<24);
        uint32_t n_sep = (uint32_t)tmp[9] | ((uint32_t)tmp[10]<<8) | ((uint32_t)tmp[11]<<16) | ((uint32_t)tmp[12]<<24);
        if((size_t)13 + n_sep + bin_len > tdl){ free(tmp); return -1; }
        size_t restore = (size_t)orig_len;
        if(*dl < restore){ free(tmp); return -1; }
        /* `restore` is a number the payload declares; the loop below writes
         * 2*bin_len characters plus one per separator-table entry, and only
         * the final `o != restore` catches a mismatch — long after the writes
         * have happened.  A payload declaring a short text length alongside a
         * long packed buffer therefore walked past the caller's capacity.
         * Bound every write instead of trusting a relation: n_sep is the
         * table's byte length, not its entry count, so no closed-form check
         * is available anyway. */
        const size_t dcap = *dl;
        const uint8_t *bin = tmp + 13;
        const uint8_t *tab = bin + bin_len;
        static const char *hexd_lc = "0123456789abcdef";
        static const char *hexd_uc = "0123456789ABCDEF";
        const char *hexd = hex_case ? hexd_uc : hexd_lc;
        size_t o = 0, oi = 0, vi = 0;

        for(size_t i = 0; i < bin_len; i++) {
            for(;;) {                        /* insert chars at value index vi */
                size_t save = oi;
                if(oi >= n_sep) break;
                uint32_t idx = 0; int sh = 0;
                for (;;) {
                    if (sh > 28) { free(tmp); return -1; }
                    uint8_t b = tab[oi++];
                    idx |= (uint32_t)(b & 0x7F) << sh;
                    if (!(b & 0x80)) break;
                    sh += 7;
                    if (oi >= n_sep && (b & 0x80)) { free(tmp); return -1; }
                }
                if(idx == vi) {
                    if (o + 1 > dcap) { free(tmp); return -1; }
                    d[o++] = tab[oi++];
                }
                else { oi = save; break; }
            }
            if (o + 2 > dcap) { free(tmp); return -1; }
            d[o++] = (uint8_t)hexd[bin[i] >> 4];
            d[o++] = (uint8_t)hexd[bin[i] & 0xF];
            vi++;
        }
        while(oi < n_sep) {                  /* trailing entries (idx == bin_len) */
            uint32_t idx = 0; int sh = 0;
            for (;;) {
                if (oi >= n_sep || sh > 28) { free(tmp); return -1; }
                uint8_t b = tab[oi++];
                idx |= (uint32_t)(b & 0x7F) << sh;
                if (!(b & 0x80)) break;
                sh += 7;
            }
            if(idx != vi){ free(tmp); return -1; }
            if (o + 1 > dcap) { free(tmp); return -1; }
            d[o++] = tab[oi++];
        }
        if(o != restore){ free(tmp); return -1; }
        *dl = o;
        free(tmp);
        return 0;
    }
    if(h->stride==CPENC_STRIDE_BWT){
        if(sl < CPENC_HEADER_SIZE+4) return -1;
        const uint8_t *rp = s + CPENC_HEADER_SIZE;
        uint32_t oi=(uint32_t)rp[0]|((uint32_t)rp[1]<<8)|((uint32_t)rp[2]<<16)|((uint32_t)rp[3]<<24);
        size_t pl = h->compressed_size - 4;
        if(pl > (size_t)(sl - CPENC_HEADER_SIZE - 4)) return -1;
        /* Redundant with the entry check above, and kept on purpose: this is
         * the branch that lost it once, and the memcpy at the end writes
         * exactly this many bytes. */
        if(*dl < h->original_size) return -1;
        uint8_t *bw = malloc(h->original_size + 1);
        if(!bw) return -2;
        size_t outl = h->original_size + 1;
        int rc = cpenc_decompress(rp + 4, pl, bw, &outl);
        if(rc != 0 || outl != h->original_size + 1){ free(bw); return -1; }
        mtf_dec_inplace(bw, outl);               /* undo MTF -> BWT stream */
        uint8_t *tmp = malloc(outl);
        if(!tmp){ free(bw); return -2; }
        /* n+1 rows, the EOF row first.  bwt_dec validates the primary index
         * and the LF chain: both come out of the untrusted payload. */
        if (bwt_dec(bw, outl, oi, tmp) != 0) { free(tmp); free(bw); return -1; }
        memcpy(d, tmp + 1, h->original_size);    /* strip the EOF symbol */
        free(tmp);
        free(bw);
        *dl = h->original_size;
        return 0;
    }

    if(h->stride==CPENC_STRIDE_BWT2){
        /* v7: [u8 k][eof_row u32][k * u32 anchors][inner cpenc block].
         * The inner block decodes to the RLE2 stream; undo RLE2, undo MTF,
         * then walk the LF cycle as k interleaved chains and strip the
         * sentinel exactly as the 0xFC path does.  Everything in the
         * payload is untrusted: k, the primary index, the anchors and the
         * inner header all come off disk. */
        if(sl < CPENC_HEADER_SIZE + 1 + 4) return -1;
        const uint8_t *rp = s + CPENC_HEADER_SIZE;
        size_t pl = h->compressed_size;
        if(pl > (size_t)(sl - CPENC_HEADER_SIZE)) return -1;
        unsigned k = rp[0];
        if(k < 8 || k > 16) return -1;           /* the format's own range */
        size_t ovh = 1 + 4 + (size_t)k * 4;
        if(pl < ovh + CPENC_HEADER_SIZE) return -1;
        /* Redundant with the entry check, as in the 0xFC branch: this is
         * the memcpy that writes original_size bytes. */
        if(*dl < h->original_size) return -1;
        size_t n1 = (size_t)h->original_size + 1;
        if(n1 < (size_t)k) return -1;
        uint32_t oi = (uint32_t)rp[1] | ((uint32_t)rp[2] << 8) |
                      ((uint32_t)rp[3] << 16) | ((uint32_t)rp[4] << 24);
        uint32_t anc[16];
        for (unsigned j = 0; j < k; j++) {
            const uint8_t *q = rp + 5 + (size_t)j * 4;
            anc[j] = (uint32_t)q[0] | ((uint32_t)q[1] << 8) |
                     ((uint32_t)q[2] << 16) | ((uint32_t)q[3] << 24);
        }
        /* Reject a corrupt table before allocating: the walk below must
         * never start from a position it cannot trust.  The anchors are
         * positions along the LF cycle, so they are in range but unordered;
         * only duplicates are structurally invalid (two chains would start
         * on the same cycle position).  bwt_dec_anchors repeats these
         * checks (defence in depth) and bounds every step. */
        if(oi >= n1 || anc[0] != 0) return -1;
        for (unsigned j = 0; j < k; j++) {
            if(anc[j] >= n1) return -1;
            for (unsigned i = 0; i < j; i++)
                if(anc[j] == anc[i]) return -1;
        }
        const uint8_t *ip = rp + ovh;
        size_t ilen = pl - ovh;
        /* The inner block's declared size is the RLE2 stream length.  RLE2
         * output is at most 2 bytes per input byte, so a stream longer than
         * 2*n1 cannot belong to this block — and this bound keeps the
         * malloc from following a corrupt header anywhere.  HEX/SEG inner
         * blocks restore to more than they declare, so those get the
         * expanded bound; the recursive decoder checks every write against
         * the capacity it is handed. */
        uint32_t in_orig = (uint32_t)ip[4] | ((uint32_t)ip[5] << 8) |
                           ((uint32_t)ip[6] << 16) | ((uint32_t)ip[7] << 24);
        if(in_orig == 0 || (size_t)in_orig > 2 * n1) return -1;
        size_t rcap = (ip[3] == CPENC_STRIDE_HEX || ip[3] == CPENC_STRIDE_SEG)
                      ? 4 * n1 + 64 : (size_t)in_orig;
        uint8_t *r = malloc(rcap);
        uint8_t *bw = malloc(n1);
        uint8_t *tmp = malloc(n1);
        if(!r || !bw || !tmp){ free(r); free(bw); free(tmp); return -2; }
        size_t outl = rcap;
        int rc = cpenc_decompress(ip, ilen, r, &outl);
        if(rc != 0 || outl > rcap) rc = -1;
        /* Undo RLE2 and MTF in one pass (no intermediate buffer), then walk
         * the LF cycle as k interleaved chains. */
        else if(rle2_mtf_dec(r, outl, bw, n1) != n1) rc = -1;
        else rc = bwt_dec_anchors(bw, n1, oi, anc, (int)k, tmp);
        if(rc == 0){
            memcpy(d, tmp + 1, h->original_size);   /* strip the EOF symbol */
            *dl = h->original_size;
        }
        free(r); free(bw); free(tmp);
        return rc == 0 ? 0 : -1;
    }

    if(h->stride==CPENC_STRIDE_REF){
        /* payload = compressed token stream; its own header carries the
         * stream length, which is what ref_decode needs to know */
        if(sl < CPENC_HEADER_SIZE + CPENC_HEADER_SIZE) return -1;
        if(*dl < h->original_size) return -1;
        const uint8_t *rp = s + CPENC_HEADER_SIZE;
        /* The stream length comes from the same untrusted header as
         * everything else, so it has to be checked against what is actually
         * there before it is handed to the inner decoder as a buffer length:
         * a header claiming more than the block holds made the inner decode
         * read past the end (fuzzing caught it; REF2 already had this check
         * and REF did not). */
        size_t pln = h->compressed_size;
        if(pln > (size_t)(sl - CPENC_HEADER_SIZE) || pln < CPENC_HEADER_SIZE) return -1;
        uint32_t tn = (uint32_t)rp[4] | ((uint32_t)rp[5]<<8) |
                      ((uint32_t)rp[6]<<16) | ((uint32_t)rp[7]<<24);
        /* tn comes from the payload: bound it against the output it can
         * possibly produce before it becomes a malloc size. */
        if(!tn || tn > (size_t)h->original_size * 4 + 64) return -1;
        uint8_t *tok = malloc(tn);
        if(!tok) return -2;
        size_t outl = tn;
        if(cpenc_decompress(rp, pln, tok, &outl) != 0 || outl != tn){
            free(tok); return -1;
        }
        int rc = ref_decode(tok, tn, d, h->original_size);
        free(tok);
        if(rc != 0) return -1;
        *dl = h->original_size;
        return 0;
    }

    if(h->stride==CPENC_STRIDE_REF2){
        /* v2: columnar sub-streams.  Rebuild the v1 token stream and
         * hand it to the unchanged ref_decode. */
        if(sl < CPENC_HEADER_SIZE + REF2_HDR) return -1;
        if(*dl < h->original_size) return -1;
        const uint8_t *pl = s + CPENC_HEADER_SIZE;
        size_t pln = h->compressed_size;
        if(pln > (size_t)(sl - CPENC_HEADER_SIZE) || pln < REF2_HDR) return -1;
        uint32_t tn = (uint32_t)pl[0] | ((uint32_t)pl[1]<<8) |
                      ((uint32_t)pl[2]<<16) | ((uint32_t)pl[3]<<24);
        /* tn comes from the payload: bound it against the output it can
         * possibly produce before it becomes a malloc size. */
        if(!tn || tn > (size_t)h->original_size * 4 + 64) return -1;
        uint8_t *tok = malloc(tn);
        if(!tok) return -2;
        if(ref_unpack2(pl, pln, tok, tn) != 0){ free(tok); return -1; }
        int rc = ref_decode(tok, tn, d, h->original_size);
        free(tok);
        if(rc != 0) return -1;
        *dl = h->original_size;
        return 0;
    }

    if(h->stride==CPENC_STRIDE_RAW){
        if(h->original_size==0){*dl=0;return 0;}
        if(*dl<h->original_size) return -1;
        const uint8_t *rp = s+CPENC_HEADER_SIZE;
        int lag_idx = CPENC_LAG_NONE;
        if(h->version>=3){
            if(sl<CPENC_HEADER_SIZE+1+(size_t)h->original_size) return -1;
            lag_idx=*rp++;
        } else if((size_t)(sl-CPENC_HEADER_SIZE)<(size_t)h->original_size) return -1;
        memcpy(d,rp,h->original_size);
        if(lag_idx){
            for(size_t i=(size_t)lag_idx;i<h->original_size;i++) d[i]+=d[i-(size_t)lag_idx];
        }
        *dl=h->original_size;return 0;}

    int stride=(int)h->stride; if(stride<1||stride>8) return -1;
    size_t orig=h->original_size, gs=orig/(size_t)stride;
    if(*dl<orig) return -1;
    const uint8_t *p=s+CPENC_HEADER_SIZE, *edd=s+sl;
    int lag_idx = CPENC_LAG_NONE;
    if(h->version>=3){ if(p>=edd)return-1; lag_idx=*p++; }
    uint8_t orders[8];
    for(int g=0;g<stride;g++){if(p>=edd)return-1;orders[g]=*p++;}

    uint8_t *trans = (stride==1) ? d : malloc(orig);
    if(!trans) return -2;
    int fused[8] = {0};   /* groups whose un-delta the fused decoder already did */

    for(int g=0;g<stride;g++){
        if(p>=edd){if(stride>1)free(trans);return-1;}
        uint8_t tag=*p++;
        uint8_t *gp=trans+g*gs;

        if(tag==3){
            if(REMAIN(p,edd)<4){if(stride>1)free(trans);return-1;}
            uint32_t rl=(uint32_t)p[0]|((uint32_t)p[1]<<8)|((uint32_t)p[2]<<16)|((uint32_t)p[3]<<24);p+=4;
            if(rl!=gs||REMAIN(p,edd)<rl){if(stride>1)free(trans);return-1;}
            memcpy(gp,p,rl);p+=rl;
        }else if(tag==0){
            if(REMAIN(p,edd)<6){if(stride>1)free(trans);return-1;}
            uint32_t tgs=(uint32_t)p[0]|((uint32_t)p[1]<<8)|((uint32_t)p[2]<<16)|((uint32_t)p[3]<<24);p+=4;
            uint8_t val=*p++;p++;
            if(tgs!=gs){if(stride>1)free(trans);return-1;}
            memset(gp,val,gs);
        }else if(tag==1){
            if(REMAIN(p,edd)<6){if(stride>1)free(trans);return-1;}
            uint32_t tgs=(uint32_t)p[0]|((uint32_t)p[1]<<8)|((uint32_t)p[2]<<16)|((uint32_t)p[3]<<24);p+=4;
            uint8_t nz1=*p++;int nz=nz1+1;
            if(nz>256){if(stride>1)free(trans);return-1;}
            int nbits=1;while((1<<nbits)<nz&&nbits<HUFF_MAX)nbits++;
            uint8_t lut[256];
            for(int i=0;i<nz;i++){if(p>=edd){if(stride>1)free(trans);return-1;}lut[i]=*p++;}
            if(tgs!=gs){if(stride>1)free(trans);return-1;}
            /* Use fast 32-bit batch decode when nbits <= 8 */
            if(nbits <= 8){
                uint32_t mask = (uint32_t)((1 << nbits) - 1);
                uint64_t rbuf = 0; int rbc = 0;
                size_t i = 0;
#if defined(__aarch64__)
                /* NEON: 16 symbols per batch (16*nbits bits = 2*nbits whole
                 * bytes, so every batch starts byte-aligned).  Covers
                 * nz<=64, which is every LUT group in the reference set. */
                if (nz <= 64) {
                    uint8_t t64[64];
                    for (int k = 0; k < nz; k++) t64[k] = lut[k];
                    for (int k = nz; k < 64; k++) t64[k] = lut[0];
                    for (; i + 16 <= gs && REMAIN(p,edd) >= 17; i += 16) {
                        vst1q_u8(gp + i, ne_lut16(p, nbits, t64));
                        p += 2 * nbits;
                    }
                }
#endif
#if defined(__x86_64__) && defined(__AVX2__)
                /* AVX2: same 16-symbol batch contract as the NEON path, so
                 * the bounds, the 2*nbits stride and the scalar tail below
                 * are unchanged. */
                if (nz <= 64) {
                    uint8_t t64[64];
                    for (int k = 0; k < nz; k++) t64[k] = lut[k];
                    for (int k = nz; k < 64; k++) t64[k] = lut[0];
                    for (; i + 16 <= gs && REMAIN(p,edd) >= 17; i += 16) {
                        _mm_storeu_si128((__m128i *)(gp + i), avx2_lut16(p, nbits, t64));
                        p += 2 * nbits;
                    }
                }
#endif
                /* Process 4 bytes (32 bits) at a time where possible */
                for (; i + 4 <= gs; i += 4) {
                    /* Fill exactly the 4 symbols' worth of bits: the old
                     * nbits*4+8 target drifted one byte per group (each
                     * batch consumed nbits*4 but buffered nbits*4+8),
                     * corrupting every following group. */
                    while (rbc < nbits * 4 && p < edd) { rbuf |= (uint64_t)(*p++) << rbc; rbc += 8; }
                    uint8_t i0 = rbuf & mask; rbuf >>= nbits; rbc -= nbits;
                    uint8_t i1 = rbuf & mask; rbuf >>= nbits; rbc -= nbits;
                    uint8_t i2 = rbuf & mask; rbuf >>= nbits; rbc -= nbits;
                    uint8_t i3 = rbuf & mask; rbuf >>= nbits; rbc -= nbits;
                    gp[i]   = i0 < nz ? lut[i0] : lut[0];
                    gp[i+1] = i1 < nz ? lut[i1] : lut[0];
                    gp[i+2] = i2 < nz ? lut[i2] : lut[0];
                    gp[i+3] = i3 < nz ? lut[i3] : lut[0];
                }
                for (; i < gs; i++) {
                    while (rbc < nbits && p < edd) { rbuf |= (uint64_t)(*p++) << rbc; rbc += 8; }
                    uint8_t idx = rbuf & mask; rbuf >>= nbits; rbc -= nbits;
                    gp[i] = idx < nz ? lut[idx] : lut[0];
                }
            } else {
                /* nbits > 8: use original bit-by-bit loop */
                uint64_t rbuf=0;int rbc=0;
                for(size_t i=0;i<gs;i++){
                    while(rbc<nbits&&p<edd){rbuf|=(uint64_t)(*p++)<<rbc;rbc+=8;}
                    uint8_t idx=rbuf&((1<<nbits)-1);rbuf>>=nbits;rbc-=nbits;
                    gp[i]=idx<nz?lut[idx]:lut[0];
                }
            }
        }else if(tag==6){
            /* bit planes: [8 markers][per plane: 1 value byte | pn packed
             * bytes | 4-byte count + (gap LEB, byte)*] */
            if(REMAIN(p,edd)<8){if(stride>1)free(trans);return-1;}
            const uint8_t *mk=p; p+=8;
            size_t pn_bp=(gs+7)/8;
            memset(gp,0,gs);
            for(int b=0;b<8;b++){
                if(mk[b]==0){
                    if(p>=edd){if(stride>1)free(trans);return-1;}
                    if(*p++) for(size_t i2=0;i2<gs;i2++) gp[i2]|=(uint8_t)(1u<<b);
                }else if(mk[b]==1){
                    if(REMAIN(p,edd)<pn_bp){if(stride>1)free(trans);return-1;}
                    for(size_t i2=0;i2<gs;i2++)
                        if((p[i2>>3]>>(i2&7))&1) gp[i2]|=(uint8_t)(1u<<b);
                    p+=pn_bp;
                }else{
                    if(REMAIN(p,edd)<4){if(stride>1)free(trans);return-1;}
                    uint32_t nz=(uint32_t)p[0]|((uint32_t)p[1]<<8)|((uint32_t)p[2]<<16)|((uint32_t)p[3]<<24);p+=4;
                    uint32_t prev=0;
                    for(uint32_t k=0;k<nz;k++){
                        uint32_t g2=0,sh2=0;
                        for(;;){
                            if(p>=edd||sh2>=32){if(stride>1)free(trans);return-1;}
                            uint8_t bb=*p++;
                            g2|=(uint32_t)(bb&0x7F)<<sh2;
                            if(!(bb&0x80))break;
                            sh2+=7;
                        }
                        prev+=g2;
                        if(prev>=pn_bp||p>=edd){if(stride>1)free(trans);return-1;}
                        uint8_t v=*p++;
                        for(int j2=0;j2<8;j2++) if((v>>j2)&1){
                            size_t idx=(size_t)prev*8+j2;
                            if(idx<gs) gp[idx]|=(uint8_t)(1u<<b);
                        }
                    }
                    if(mk[b]==3) for(size_t i2=0;i2<gs;i2++) gp[i2]^=(uint8_t)(1u<<b);
                }
            }
        }else if(tag==5){
            /* const + exceptions: [val][nex][pos LEB128][values] */
            if(REMAIN(p,edd)<6){if(stride>1)free(trans);return-1;}
            uint16_t cval=(uint16_t)(p[0]|(p[1]<<8));p+=2;
            uint32_t nex=(uint32_t)p[0]|((uint32_t)p[1]<<8)|((uint32_t)p[2]<<16)|((uint32_t)p[3]<<24);p+=4;
            if(nex>gs){if(stride>1)free(trans);return-1;}
            memset(gp,cval,gs);
            uint32_t pos=0;
            for(uint32_t k=0;k<nex;k++){
                uint32_t gap=0,shift=0;
                for(;;){
                    if(p>=edd||shift>28){if(stride>1)free(trans);return-1;}
                    uint8_t b=*p++;
                    gap|=(uint32_t)(b&0x7F)<<shift;
                    if(!(b&0x80))break;
                    shift+=7;
                }
                pos = (k == 0) ? gap : (pos + 1 + gap);
                if(p>=edd||pos>=gs){if(stride>1)free(trans);return-1;}
                gp[pos]=*p++;
            }
        }else if(tag==4){
            /* RLE: decode runs and vals sub-blocks, alternate fill */
            if(p>=edd){if(stride>1)free(trans);return-1;}
            uint8_t st=*p++;
            size_t rb_cap = rle_max_run_bytes(gs);   /* encoded run bytes */
            uint32_t *runs_buf = malloc(rb_cap*sizeof(uint32_t));
            uint8_t *vals_buf = malloc(gs);
            uint8_t *scratch = malloc(gs);
            if(!runs_buf||!vals_buf||!scratch){free(runs_buf);free(vals_buf);free(scratch);if(stride>1)free(trans);return-1;}
            size_t nr=0,nv=0;
            if(st==0){
                if(REMAIN(p,edd)<6){free(runs_buf);free(vals_buf);free(scratch);if(stride>1)free(trans);return-1;}
                uint32_t tgs=(uint32_t)p[0]|((uint32_t)p[1]<<8)|((uint32_t)p[2]<<16)|((uint32_t)p[3]<<24);p+=4;
                uint8_t val=*p++;p++;
                if(tgs>rb_cap){free(runs_buf);free(vals_buf);free(scratch);if(stride>1)free(trans);return-1;}
                for(size_t k=0;k<tgs;k++)runs_buf[nr++]=(uint32_t)val;
            }else if(st==1){
                if(REMAIN(p,edd)<6){free(runs_buf);free(vals_buf);free(scratch);if(stride>1)free(trans);return-1;}
                uint32_t tgs=(uint32_t)p[0]|((uint32_t)p[1]<<8)|((uint32_t)p[2]<<16)|((uint32_t)p[3]<<24);p+=4;
                uint8_t nz1=*p++;int nz=nz1+1;
                uint8_t lut[256];for(int i=0;i<nz;i++){if(p>=edd){free(runs_buf);free(vals_buf);free(scratch);if(stride>1)free(trans);return-1;}lut[i]=*p++;}
                if(tgs>rb_cap){free(runs_buf);free(vals_buf);free(scratch);if(stride>1)free(trans);return-1;}
                int nbits=1;while((1<<nbits)<nz&&nbits<12)nbits++;
                uint64_t rbuf=0;int rbc=0;
                size_t i=0;
#if defined(__aarch64__)
                if (nbits <= 8 && nz <= 64) {
                    uint8_t t64[64], tmp[16];
                    for (int k = 0; k < nz; k++) t64[k] = lut[k];
                    for (int k = nz; k < 64; k++) t64[k] = lut[0];
                    for (; i + 16 <= tgs && REMAIN(p,edd) >= 17; i += 16) {
                        vst1q_u8(tmp, ne_lut16(p, nbits, t64));
                        for (int j = 0; j < 16; j++) runs_buf[nr++] = (uint32_t)tmp[j];
                        p += 2 * nbits;
                    }
                }
#endif
#if defined(__x86_64__) && defined(__AVX2__)
                if (nbits <= 8 && nz <= 64) {
                    uint8_t t64[64], tmp[16];
                    for (int k = 0; k < nz; k++) t64[k] = lut[k];
                    for (int k = nz; k < 64; k++) t64[k] = lut[0];
                    for (; i + 16 <= tgs && REMAIN(p,edd) >= 17; i += 16) {
                        _mm_storeu_si128((__m128i *)tmp, avx2_lut16(p, nbits, t64));
                        for (int j = 0; j < 16; j++) runs_buf[nr++] = (uint32_t)tmp[j];
                        p += 2 * nbits;
                    }
                }
#endif
                for(;i<tgs;i++){
                    while(rbc<nbits&&p<edd){rbuf|=(uint64_t)(*p++)<<rbc;rbc+=8;}
                    uint8_t idx=rbuf&((1<<nbits)-1);rbuf>>=nbits;rbc-=nbits;
                    runs_buf[nr++]=(uint32_t)(idx<nz?lut[idx]:lut[0]);
                }
            }else if(st==2){
                /* huffman sub-block: symbol count lives in the stream header */
                if(REMAIN(p,edd)<264){free(runs_buf);free(vals_buf);free(scratch);if(stride>1)free(trans);return-1;}
                uint32_t rn=(uint32_t)p[256]|((uint32_t)p[257]<<8)|((uint32_t)p[258]<<16)|((uint32_t)p[259]<<24);
                /* Two different bounds: rn counts entries in runs_buf (sized
                 * rb_cap) and bytes written into scratch (sized gs).  For a
                 * small group rb_cap is the larger of the two — its +8 slack
                 * dominates — so checking only rb_cap let h_dec write past
                 * scratch. */
                if(rn>rb_cap||rn>gs){free(runs_buf);free(vals_buf);free(scratch);if(stride>1)free(trans);return-1;}
                if(h_dec(&p,edd,scratch,rn)<0){free(runs_buf);free(vals_buf);free(scratch);if(stride>1)free(trans);return-1;}
                for(size_t k=0;k<rn;k++)runs_buf[nr++]=(uint32_t)scratch[k];
            }else if(st==3){
                if(REMAIN(p,edd)<4){free(runs_buf);free(vals_buf);free(scratch);if(stride>1)free(trans);return-1;}
                uint32_t rl=(uint32_t)p[0]|((uint32_t)p[1]<<8)|((uint32_t)p[2]<<16)|((uint32_t)p[3]<<24);p+=4;
                if(rl>rb_cap||REMAIN(p,edd)<rl){free(runs_buf);free(vals_buf);free(scratch);if(stride>1)free(trans);return-1;}
                for(size_t k=0;k<rl;k++)runs_buf[nr++]=(uint32_t)p[k];
                p+=rl;
            }else if(st==4){
                uint32_t tlen=0;
                /* Capacity is scratch's gs, not rb_cap: the two differ for
                 * small groups, and ASan caught the memset in sub5_decode
                 * writing two bytes past a 4-byte scratch. */
                if(sub5_decode(&p,edd,scratch,gs,&tlen)!=0){free(runs_buf);free(vals_buf);free(scratch);if(stride>1)free(trans);return-1;}
                if(tlen>rb_cap){free(runs_buf);free(vals_buf);free(scratch);if(stride>1)free(trans);return-1;}
                for(size_t k=0;k<tlen;k++)runs_buf[nr++]=(uint32_t)scratch[k];
            }
            /* vals sub-block */

            if(p>=edd){free(runs_buf);free(vals_buf);free(scratch);if(stride>1)free(trans);return-1;}
            st=*p++;

            if(st==0){
                if(REMAIN(p,edd)<6){free(runs_buf);free(vals_buf);free(scratch);if(stride>1)free(trans);return-1;}
                uint32_t tgs=(uint32_t)p[0]|((uint32_t)p[1]<<8)|((uint32_t)p[2]<<16)|((uint32_t)p[3]<<24);p+=4;
                uint8_t val=*p++;p++;
                if(tgs>gs){free(runs_buf);free(vals_buf);free(scratch);if(stride>1)free(trans);return-1;}
                for(size_t k=0;k<tgs;k++)vals_buf[nv++]=(uint8_t)val;
            }else if(st==1){
                if(REMAIN(p,edd)<6){free(runs_buf);free(vals_buf);free(scratch);if(stride>1)free(trans);return-1;}
                uint32_t tgs=(uint32_t)p[0]|((uint32_t)p[1]<<8)|((uint32_t)p[2]<<16)|((uint32_t)p[3]<<24);p+=4;
                uint8_t nz1=*p++;int nz=nz1+1;
                uint8_t lut[256];for(int i=0;i<nz;i++){if(p>=edd){free(runs_buf);free(vals_buf);free(scratch);if(stride>1)free(trans);return-1;}lut[i]=*p++;}
                if(tgs>gs){free(runs_buf);free(vals_buf);free(scratch);if(stride>1)free(trans);return-1;}
                int nbits=1;while((1<<nbits)<nz&&nbits<12)nbits++;
                uint64_t rbuf=0;int rbc=0;
                size_t i=0;
#if defined(__aarch64__)
                if (nbits <= 8 && nz <= 64) {
                    uint8_t t64[64];
                    for (int k = 0; k < nz; k++) t64[k] = lut[k];
                    for (int k = nz; k < 64; k++) t64[k] = lut[0];
                    for (; i + 16 <= tgs && REMAIN(p,edd) >= 17; i += 16) {
                        vst1q_u8(vals_buf + nv, ne_lut16(p, nbits, t64));
                        nv += 16;
                        p += 2 * nbits;
                    }
                }
#endif
#if defined(__x86_64__) && defined(__AVX2__)
                if (nbits <= 8 && nz <= 64) {
                    uint8_t t64[64];
                    for (int k = 0; k < nz; k++) t64[k] = lut[k];
                    for (int k = nz; k < 64; k++) t64[k] = lut[0];
                    for (; i + 16 <= tgs && REMAIN(p,edd) >= 17; i += 16) {
                        _mm_storeu_si128((__m128i *)(vals_buf + nv), avx2_lut16(p, nbits, t64));
                        nv += 16;
                        p += 2 * nbits;
                    }
                }
#endif
                for(;i<tgs;i++){
                    while(rbc<nbits&&p<edd){rbuf|=(uint64_t)(*p++)<<rbc;rbc+=8;}
                    uint8_t idx=rbuf&((1<<nbits)-1);rbuf>>=nbits;rbc-=nbits;
                    vals_buf[nv++]=(uint8_t)(idx<nz?lut[idx]:lut[0]);
                }
            }else if(st==2){
                if(REMAIN(p,edd)<264){free(runs_buf);free(vals_buf);free(scratch);if(stride>1)free(trans);return-1;}
                uint32_t rn=(uint32_t)p[256]|((uint32_t)p[257]<<8)|((uint32_t)p[258]<<16)|((uint32_t)p[259]<<24);
                if(rn>gs){free(runs_buf);free(vals_buf);free(scratch);if(stride>1)free(trans);return-1;}
                if(h_dec(&p,edd,vals_buf,rn)<0){free(runs_buf);free(vals_buf);free(scratch);if(stride>1)free(trans);return-1;}
                nv=rn;
            }else if(st==3){
                if(REMAIN(p,edd)<4){free(runs_buf);free(vals_buf);free(scratch);if(stride>1)free(trans);return-1;}
                uint32_t rl=(uint32_t)p[0]|((uint32_t)p[1]<<8)|((uint32_t)p[2]<<16)|((uint32_t)p[3]<<24);p+=4;
                if(rl>gs||REMAIN(p,edd)<rl){free(runs_buf);free(vals_buf);free(scratch);if(stride>1)free(trans);return-1;}
                memcpy(vals_buf,p,rl);nv=rl;p+=rl;
            }else if(st==4){
                uint32_t tlen=0;
                if(sub5_decode(&p,edd,vals_buf,gs,&tlen)!=0){free(runs_buf);free(vals_buf);free(scratch);if(stride>1)free(trans);return-1;}
                nv=tlen;
            }
            /* alternate fill: run zeros, then value */
            size_t di=0,ri=0,vi=0;
            while(di<gs){
                /* decode one LEB128 run length */
                uint32_t r=0,shift=0;

                for(;;){
                    if(ri>=nr){free(runs_buf);free(vals_buf);free(scratch);if(stride>1)free(trans);return-1;}
                    uint8_t b=(uint8_t)runs_buf[ri++];
                    if(b<255){r=(uint32_t)b;break;}
                    r=255;
                    /* escape: remainder in LEB128 — accumulate with + (r is
                     * preset 255 whose bits are all set; |= would no-op) */
                    for(;;){
                        if(ri>=nr||shift>28){free(runs_buf);free(vals_buf);free(scratch);if(stride>1)free(trans);return-1;}
                        uint8_t c2=(uint8_t)runs_buf[ri++];
                        r+=(uint32_t)(c2&0x7F)<<shift;
                        if(!(c2&0x80))break;
                        shift+=7;
                    }
                    break;
                }
                if(r>gs-di)r=(uint32_t)(gs-di);
                memset(gp+di,0,r);di+=r;
                if(di<gs){
                    if(vi>=nv){free(runs_buf);free(vals_buf);free(scratch);if(stride>1)free(trans);return-1;}
                    gp[di++]=vals_buf[vi++];
                }
            }
            free(runs_buf);free(vals_buf);free(scratch);
        }else{
            int hrc;
            if(stride==1 && orders[g]==1 && lag_idx==0){
                hrc = h_dec_p(&p,edd,gp,gs,1);       /* fused prefix sum */
                fused[g] = 1;
            }
            else if(gs >= (1u<<20)) hrc = h_dec_seg(&p,edd,gp,gs);
            else                     hrc = h_dec(&p,edd,gp,gs);
            if(hrc<0){if(stride>1)free(trans);return-1;}
        }
    }
    /* v2: read the stride remainder tail bytes (never delta-coded) */
    if (h->version >= 2) {
        size_t gend = gs * (size_t)stride;
        if (orig > gend) {
            size_t rem = orig - gend;
            if (REMAIN(p,edd) < rem) { if (stride > 1) free(trans); return -1; }
            memcpy(trans + gend, p, rem);
            p += rem;
        }
    }
    /* Un-delta each group.  A group decoded through the fused prefix-sum
     * decoder (h_dec_p) already has the un-delta folded into its output,
     * so un-deltaing it here a second time corrupts the block — the
     * stride=1/order=1/Huffman case.  h_dec_p is only reached when
     * stride==1, so the fused group is group 0. */
    for(int g=0;g<stride;g++)
        if(!fused[g]) dd(trans+g*gs,gs,orders[g]);
    if(stride>1){ur(trans,d,orig,stride);free(trans);}
    /* v3: undo lag-delta (applied pre-transpose on the raw byte order) */
    if(lag_idx){
        for(size_t i=(size_t)lag_idx;i<orig;i++) d[i]+=d[i-(size_t)lag_idx];
    }
    *dl=orig; return 0;
}

int cpenc_compress(const void *sv, size_t sl, void *dv, size_t *dl, int mo) {
    return cpenc_compress_impl(sv, sl, dv, dl, mo, 0, 0);
}
