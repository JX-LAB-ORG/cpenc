#include "bwt.h"
#include <stdlib.h>
#include <string.h>

typedef struct { uint32_t idx, r1, r2; } Sfx;

/* LSD radix on one 32-bit key field: low 16 bits first, then high 16
 * (high-first would be MSD, which needs bucket recursion — wrong here). */
static void radix_pass(Sfx *a, Sfx *tmp, size_t n, int hi16) {
    /* LSD order: low 16 bits first (stable), then high 16 bits. */
    uint32_t cnt[65536];
    memset(cnt, 0, sizeof(cnt));
    for (size_t i = 0; i < n; i++) {
        uint32_t v = hi16 ? (a[i].r1 & 0xFFFF) : (a[i].r2 & 0xFFFF);
        cnt[v]++;
    }
    uint32_t sum = 0;
    for (int i = 0; i < 65536; i++) { uint32_t c = cnt[i]; cnt[i] = sum; sum += c; }
    for (size_t i = 0; i < n; i++) {
        uint32_t v = hi16 ? (a[i].r1 & 0xFFFF) : (a[i].r2 & 0xFFFF);
        tmp[cnt[v]++] = a[i];
    }
    memset(cnt, 0, sizeof(cnt));
    for (size_t i = 0; i < n; i++) {
        uint32_t v = hi16 ? (tmp[i].r1 >> 16) : (tmp[i].r2 >> 16);
        cnt[v]++;
    }
    sum = 0;
    for (int i = 0; i < 65536; i++) { uint32_t c = cnt[i]; cnt[i] = sum; sum += c; }
    for (size_t i = 0; i < n; i++) {
        uint32_t v = hi16 ? (tmp[i].r1 >> 16) : (tmp[i].r2 >> 16);
        a[cnt[v]++] = tmp[i];
    }
}

size_t bwt_enc(const uint8_t *s, size_t n, uint8_t *out) {
    if (n == 0) return 0;
    uint32_t *rank = malloc(n * sizeof(uint32_t));
    uint32_t *nrank = malloc(n * sizeof(uint32_t));
    Sfx *sfx = malloc(n * sizeof(Sfx));
    Sfx *tmp = malloc(n * sizeof(Sfx));
    for (size_t i = 0; i < n; i++) rank[i] = s[i];
    for (size_t k = 1; k < n; k <<= 1) {
        for (size_t i = 0; i < n; i++) {
            sfx[i].idx = (uint32_t)i;
            /* cyclic BWT (no EOF sentinel): n rows, LF-mapping decode is
             * closed on the n-row table — the decoder depends on this */
            sfx[i].r1 = rank[i];
            sfx[i].r2 = rank[(i + k) % n];
        }
        radix_pass(sfx, tmp, n, 0);   /* LSB: r2 */
        radix_pass(sfx, tmp, n, 1);   /* MSB: r1 */
        uint32_t nr = 0;
        for (size_t i = 0; i < n; i++) {
            if (i > 0 && (sfx[i].r1 != sfx[i-1].r1 || sfx[i].r2 != sfx[i-1].r2)) nr++;
            nrank[sfx[i].idx] = nr;
        }
        memcpy(rank, nrank, n * sizeof(uint32_t));
        if (nr == n - 1) break;
    }
    size_t orig = 0;
    for (size_t i = 0; i < n; i++) {
        if (sfx[i].idx == 0) orig = i;
        out[i] = s[(sfx[i].idx + n - 1) % n];
    }
    free(rank); free(nrank); free(sfx); free(tmp);
    return orig;
}

void bwt_dec(const uint8_t *b, size_t n, size_t orig_idx, uint8_t *out) {
    /* precompute the LF table in one sequential pass (start[c]+occ[pos]),
     * then the backward walk is a single table lookup per step instead of
     * three dependent loads — LF dominates BWT decode time. */
    uint32_t *cnt = calloc(256, sizeof(uint32_t));
    for (size_t i = 0; i < n; i++) cnt[b[i]]++;
    uint32_t start[256];
    uint32_t sum = 0;
    for (int i = 0; i < 256; i++) { start[i] = sum; sum += cnt[i]; }
    uint32_t *lf = malloc(n * sizeof(uint32_t));
    uint32_t cpos[256];
    memcpy(cpos, start, sizeof(cpos));
    for (size_t i = 0; i < n; i++) lf[i] = cpos[b[i]]++;
    size_t pos = orig_idx;
    for (size_t i = n; i-- > 0;) {
        out[i] = b[pos];
        pos = lf[pos];
    }
    free(cnt); free(lf);
}

void mtf_enc_inplace(uint8_t *p, size_t n) {
    uint8_t st[256];
    for (int i = 0; i < 256; i++) st[i] = (uint8_t)i;
    for (size_t i = 0; i < n; i++) {
        uint8_t v = p[i];
        int r = 0;
        while (st[r] != v) r++;
        p[i] = (uint8_t)r;
        for (int j = r; j > 0; j--) st[j] = st[j-1];
        st[0] = v;
    }
}
void mtf_dec_inplace(uint8_t *p, size_t n) {
    uint8_t st[256];
    for (int i = 0; i < 256; i++) st[i] = (uint8_t)i;
    for (size_t i = 0; i < n; i++) {
        uint8_t r = p[i];
        uint8_t v = st[r];
        p[i] = v;
        for (int j = r; j > 0; j--) st[j] = st[j-1];
        st[0] = v;
    }
}
