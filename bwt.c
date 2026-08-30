#include "bwt.h"
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

/* ── SA-IS: linear-time suffix array (Nong, Zhang, Chan 2009) ─────
 * s[0..n-1] alphabet 0..k-1, implicit sentinel s[n]=0 (strictly minimal).
 * Output sa[0..n] (n+1 suffixes; sa[0]=n is the sentinel suffix).
 * Fallbacks: strictly-decreasing input → reverse order; non-decreasing
 * recursion levels (equal chars) → insertion sort (small n only). */
static void sais_impl(const int32_t *s, int32_t *sa, int n, int k) {
    if (n <= 0) { sa[0] = 0; return; }
    uint8_t *t = malloc((size_t)n + 1);
    int32_t *bkt = malloc((size_t)(k + 1) * sizeof(int32_t));
    int32_t *cnt = malloc((size_t)k * sizeof(int32_t));
    t[n] = 1; t[n-1] = 0;
    for (int i = n-2; i >= 0; i--)
        t[i] = (s[i] < s[i+1] || (s[i] == s[i+1] && t[i+1])) ? 1 : 0;
    memset(cnt, 0, (size_t)k * sizeof(int32_t));
    for (int i = 0; i <= n; i++) cnt[s[i]]++;
    bkt[0] = 0;
    for (int i = 0; i < k; i++) bkt[i+1] = bkt[i] + cnt[i];
    int nlms = 0;
    int32_t *lms = malloc((size_t)n * sizeof(int32_t));
    for (int i = 0; i < n; i++)
        if (t[i] == 1 && (i == 0 || t[i-1] == 0)) lms[nlms++] = i;
    if (nlms == 0) {
        /* no LMS => monotone input: non-decreasing suffixes sort by
         * position (each is a prefix of the next or smaller); otherwise
         * strictly decreasing sorts reversed.  O(n), no insertion sort —
         * the O(n^2) fallback exploded on large recursion levels. */
        int non_dec = 1;
        for (int i = 0; i + 1 < n; i++)
            if (s[i] > s[i+1]) { non_dec = 0; break; }
        if (non_dec) {
            for (int i = 0; i < n; i++) sa[i+1] = (int32_t)i;
        } else {
            for (int i = 0; i < n; i++) sa[i+1] = n-1-i;
        }
        sa[0] = n;
        free(t); free(bkt); free(cnt); free(lms);
        return;
    }
    int32_t *cur = malloc((size_t)(k + 1) * sizeof(int32_t));
    for (int i = 0; i <= n; i++) sa[i] = -1;
    sa[0] = n;
    for (int i = 0; i < k; i++) cur[i] = bkt[i+1];
    for (int i = 0; i < nlms; i++) sa[--cur[s[lms[i]]]] = lms[i];
    memcpy(cur, bkt, (size_t)(k+1) * sizeof(int32_t));
    for (int i = 0; i <= n; i++) {
        int32_t p = sa[i];
        if (p > 0 && t[p-1] == 0) sa[cur[s[p-1]]++] = p-1;
    }
    for (int i = 0; i < k; i++) cur[i] = bkt[i+1];
    for (int i = n; i >= 0; i--) {
        int32_t p = sa[i];
        if (p > 0 && t[p-1] == 1) sa[--cur[s[p-1]]] = p-1;
    }
    int32_t *lms_ord = malloc((size_t)nlms * sizeof(int32_t));
    int c = 0;
    for (int i = 0; i <= n; i++)
        if (sa[i] >= 0 && sa[i] < n && t[sa[i]] == 1 &&
            (sa[i] == 0 || t[sa[i]-1] == 0)) lms_ord[c++] = sa[i];
    int32_t *name = malloc((size_t)n * sizeof(int32_t));
    int names = 0;
    int32_t prev = -1;
    for (int i = 0; i < nlms; i++) {
        int32_t p = lms_ord[i];
        if (prev >= 0) {
            int d = 0, eq = 1;
            while (1) {
                int isl_p = (d > 0) && (t[p+d] == 1 && (p+d == 0 || t[p+d-1] == 0));
                int isl_q = (d > 0) && (t[prev+d] == 1 && (prev+d == 0 || t[prev+d-1] == 0));
                if (isl_p && isl_q) break;
                if (isl_p || isl_q) { eq = 0; break; }
                if (s[p+d] != s[prev+d]) { eq = 0; break; }
                if (t[p+d] != t[prev+d]) { eq = 0; break; }
                d++;
            }
            if (!eq) names++;
        }
        name[p] = names;
        prev = p;
    }
    int32_t *s2 = malloc((size_t)(nlms + 1) * sizeof(int32_t));
    for (int i = 0; i < nlms; i++) s2[i] = name[lms[i]] + 1;
    s2[nlms] = 0;
    int32_t *sa2 = malloc((size_t)(nlms + 1) * sizeof(int32_t));
    if (names + 1 < nlms) {
        sais_impl(s2, sa2, nlms, names + 2);
    } else {
        for (int i = 0; i < nlms; i++) sa2[name[lms[i]] + 1] = i;
        sa2[0] = nlms;
    }
    for (int i = 0; i <= n; i++) sa[i] = -1;
    sa[0] = n;
    for (int i = 0; i < k; i++) cur[i] = bkt[i+1];
    for (int i = nlms; i >= 1; i--) {
        int32_t p = lms[sa2[i]];
        sa[--cur[s[p]]] = p;
    }
    memcpy(cur, bkt, (size_t)(k+1) * sizeof(int32_t));
    for (int i = 0; i <= n; i++) {
        int32_t p = sa[i];
        if (p > 0 && t[p-1] == 0) sa[cur[s[p-1]]++] = p-1;
    }
    for (int i = 0; i < k; i++) cur[i] = bkt[i+1];
    for (int i = n; i >= 0; i--) {
        int32_t p = sa[i];
        if (p > 0 && t[p-1] == 1) sa[--cur[s[p-1]]] = p-1;
    }
    free(t); free(bkt); free(cnt); free(cur);
    free(lms); free(lms_ord); free(name); free(s2); free(sa2);
}
void sais_sa(const uint8_t *s, int32_t *sa, int n) {
    int32_t *os = malloc(((size_t)n + 1) * sizeof(int32_t));
    for (int i = 0; i < n; i++) os[i] = (int32_t)s[i] + 1;
    os[n] = 0;
    sais_impl(os, sa, n, 258);
    free(os);
}

/* ── MTF (move-to-front), bzip2-style: BWT clustering → zero runs ── */
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

/* ── BWT ─────────────────────────────────────────────────────────
 * Linear BWT: suffix array of s (implicit EOF sentinel, from SA-IS).
 * The L column has n+1 rows: out[i] = predecessor of suffix sa[i],
 * where the row sa[i]==0 (the original string's start) carries the
 * EOF symbol 0.  out has n+1 entries; decode emits n+1 (last is EOF)
 * and the caller strips the 0.  Sorting half the length of the cyclic
 * s+s trick — the previous scheme ran SA-IS on 2n bytes. */
size_t bwt_enc(const uint8_t *s, size_t n, uint8_t *out) {
    if (n == 0) { out[0] = 0; return 0; }
    int32_t *sa = malloc((n + 1) * sizeof(int32_t));
    sais_sa(s, sa, n);
    for (int32_t i = 0; i <= (int32_t)n; i++) {
        int32_t p = sa[i];
        out[i] = (p == 0) ? 0 : s[p - 1];
    }
    free(sa);
    return 0;   /* EOF row is always sa[0] */
}

/* LF-mapping decode of the n+1-row linear BWT; output is n+1 bytes
 * with the EOF symbol at the last emitted position. */
void bwt_dec(const uint8_t *b, size_t n1, size_t orig_idx, uint8_t *out) {
    uint32_t *cnt = calloc(256, sizeof(uint32_t));
    for (size_t i = 0; i < n1; i++) cnt[b[i]]++;
    uint32_t start[256];
    uint32_t sum = 0;
    for (int i = 0; i < 256; i++) { start[i] = sum; sum += cnt[i]; }
    uint32_t *lf = malloc(n1 * sizeof(uint32_t));
    uint32_t cpos[256];
    memcpy(cpos, start, sizeof(cpos));
    for (size_t i = 0; i < n1; i++) lf[i] = cpos[b[i]]++;
    size_t pos = orig_idx;                    /* EOF row: sa[0] */
    for (size_t i = n1; i-- > 0;) {
        out[i] = b[pos];
        pos = lf[pos];
    }
    free(cnt); free(lf);
}
