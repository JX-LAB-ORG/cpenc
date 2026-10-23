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
        /* no LMS => monotone input: strictly increasing suffixes sort by
         * position; otherwise the order is reversed.  The test must be
         * STRICT: with equal neighbours the shorter suffix is smaller
         * (it is a prefix), so [2,2,2] sorts reversed, not by position —
         * a non-decreasing test put runs of equal bytes in the wrong
         * order and made bwt_dec output garbage. */
        int strict_inc = 1;
        for (int i = 0; i + 1 < n; i++)
            if (s[i] >= s[i+1]) { strict_inc = 0; break; }
        if (strict_inc) {
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
    size_t eof_row = 0;
    for (int32_t i = 0; i <= (int32_t)n; i++) {
        int32_t p = sa[i];
        out[i] = (p == 0) ? 0 : s[p - 1];
        if (p == 0) eof_row = (size_t)i;
    }
    free(sa);
    return eof_row;
}

/* LF-mapping decode of the n+1-row linear BWT; output is n+1 bytes with
 * the EOF symbol first.
 *
 * The EOF row is written as byte 0, which is also a legal input byte, so
 * the decoder cannot tell the sentinel apart from a real 0 in L.  It does
 * not have to: `eof_row` (returned by bwt_enc) marks it, and removing that
 * one occurrence from the character counts restores the true alphabet,
 * where the sentinel is unique and strictly smallest.  For inputs with no
 * 0x00 byte this yields exactly the same LF as the naive build. */
int bwt_dec(const uint8_t *b, size_t n1, size_t eof_row, uint8_t *out) {
    /* Both inputs come from an untrusted container.  The primary index was
     * used unchecked, and the LF chain it produces was followed unchecked:
     * a crafted block (primary index 0xffff1509 in a 16 KB block was the
     * fuzzer's example) walked off the end of b[]. */
    if (n1 == 0 || eof_row >= n1) return -1;
    uint32_t cnt[256] = {0};
    for (size_t i = 0; i < n1; i++) cnt[b[i]]++;
    cnt[b[eof_row]]--;                        /* sentinel is not a character */
    uint32_t start[256];
    uint32_t sum = 0;
    for (int i = 0; i < 256; i++) { start[i] = sum; sum += cnt[i]; }
    uint32_t *lf = malloc(n1 * sizeof(uint32_t));
    if (!lf) return -1;
    uint32_t cpos[256];
    memcpy(cpos, start, sizeof(cpos));
    for (size_t i = 0; i < n1; i++) {
        if (i == eof_row) { lf[i] = 0; continue; }   /* F[0] = sentinel */
        lf[i] = 1 + cpos[b[i]]++;
    }
    size_t pos = 0;
    for (size_t i = n1; i-- > 0;) {
        if (pos >= n1) { free(lf); return -1; }      /* chain left the block */
        out[i] = b[pos];
        pos = lf[pos];
    }
    free(lf);
    return 0;
}

/* ── anchored LF walk (v7 BWT2) ─────────────────────────────────
 * Serial decode is one dependent load per output byte.  The same LF
 * permutation is a cycle; start k chains at stored positions and the
 * loads of the k chains are independent, so one core can have several
 * cache misses in flight.
 *
 * The LF build below is deliberately the same construction as bwt_dec's:
 * cnt[] minus the EOF occurrence, start[] as the F column, lf[i] = the F
 * position of the i-th L character (the EOF row maps to F[0]).  It is
 * repeated here rather than factored out so the v6 (stride 0xFC) path
 * stays byte-for-byte the behaviour it shipped with.
 */
static uint32_t *bwt_lf_build(const uint8_t *b, size_t n1, size_t eof_row) {
    if (n1 == 0 || eof_row >= n1) return NULL;
    /* The histogram is four interleaved counters: `cnt[b[i]]++` has a
     * loop-carried dependency through the counter it touches, and four
     * independent tables (4KB of stack, L1-resident) let the core issue
     * four increments per cycle instead.  Measured 2.88 -> 2.03 ms/MB,
     * which matters because every BWT2 decode pays this. */
    uint32_t cnt[4][256];
    memset(cnt, 0, sizeof cnt);
    size_t i = 0;
    for (; i + 4 <= n1; i += 4) {
        cnt[0][b[i]]++; cnt[1][b[i + 1]]++;
        cnt[2][b[i + 2]]++; cnt[3][b[i + 3]]++;
    }
    for (; i < n1; i++) cnt[0][b[i]]++;
    for (int c = 1; c < 4; c++)
        for (int k = 0; k < 256; k++) cnt[0][k] += cnt[c][k];
    cnt[0][b[eof_row]]--;                     /* sentinel is not a character */
    uint32_t start[256];
    uint32_t sum = 0;
    for (int i = 0; i < 256; i++) { start[i] = sum; sum += cnt[0][i]; }
    uint32_t *lf = malloc(n1 * sizeof(uint32_t));
    if (!lf) return NULL;
    uint32_t cpos[256];
    memcpy(cpos, start, sizeof(cpos));
    /* Range-check each value where it is produced: this is the proof the
     * chain walk relies on, and a branch here is independent of the walk's
     * dependent loads.  (The values cannot actually exceed n1-1 for a table
     * this function built — the check exists so that a future change here
     * fails closed instead of walking outside the block.) */
    for (i = 0; i < eof_row; i++) {
        lf[i] = 1 + cpos[b[i]]++;
        if (lf[i] >= n1) { free(lf); return NULL; }
    }
    lf[eof_row] = 0;                          /* F[0] = sentinel */
    for (i = eof_row + 1; i < n1; i++) {
        lf[i] = 1 + cpos[b[i]]++;
        if (lf[i] >= n1) { free(lf); return NULL; }
    }
    return lf;
}

int bwt_make_anchors(const uint8_t *bwt, size_t n1, size_t orig_idx,
                     uint32_t *anchors, int k) {
    if (k < 1 || n1 == 0 || orig_idx >= n1) return -1;
    uint32_t *lf = bwt_lf_build(bwt, n1, orig_idx);
    if (!lf) return -1;
    size_t m = n1 / (size_t)k;
    size_t pos = 0;
    int rc = 0;
    for (int j = 0; j < k; j++) {
        anchors[j] = (uint32_t)pos;
        for (size_t t = 0; t < m; t++) {
            pos = lf[pos];
            /* A valid BWT cannot leave the block (the LF values are F
             * positions), but the walk is cheap to guard. */
            if (pos >= n1) { rc = -1; goto done; }
        }
    }
done:
    free(lf);
    return rc;
}

int bwt_dec_anchors(const uint8_t *bwt, size_t n1, size_t orig_idx,
                    const uint32_t *anchors, int k, uint8_t *out) {
    if (n1 == 0 || orig_idx >= n1) return -1;
    if (k < 8 || k > 16) return -1;              /* the v7 format's range */
    if (n1 < (size_t)k) return -1;               /* k distinct positions */
    /* The cycle starts at 0.  The other anchors are j*m LF steps along it,
     * so they are in-range but *not* ordered: the LF permutation jumps
     * around the block.  What must hold is that they are pairwise distinct;
     * a duplicate would make two chains start on the same cycle position,
     * which is a corrupt table even though each write stays in its own
     * output range. */
    if (anchors[0] != 0) return -1;
    for (int j = 0; j < k; j++) {
        if (anchors[j] >= n1) return -1;
        for (int i = 0; i < j; i++)
            if (anchors[j] == anchors[i]) return -1;
    }
    uint32_t *lf = bwt_lf_build(bwt, n1, orig_idx);
    if (!lf) return -1;
    /* bwt_lf_build proved every value in range as it wrote it, so no step
     * below can leave the block and the walk needs no per-step check: that
     * check sat on the dependency chain (each chain's next load would wait
     * for the previous comparison) and halved the walk's speed.  Together
     * with the anchor validation above, every step of every chain is in
     * bounds. */
    size_t m = n1 / (size_t)k;
    size_t rem = n1 - m * (size_t)k;             /* the last chain takes it */
    uint32_t a[16];
    uint8_t *op[16];
    memcpy(a, anchors, (size_t)k * sizeof(uint32_t));
    for (int j = 0; j < k; j++)
        op[j] = out + (n1 - 1 - (size_t)j * m);  /* each chain runs downward */
    for (size_t t = 0; t < m; t++) {
        for (int j = 0; j < k; j++) {
            uint32_t ps = a[j];
            *op[j]-- = bwt[ps];
            a[j] = lf[ps];
        }
    }
    {
        /* Chain k-1's anchored start covered steps [0,m); op[k-1] now points
         * at the byte before that range, so it continues through the
         * remainder and every output byte is written once. */
        uint8_t *q = op[k - 1];
        uint32_t ps = a[k - 1];
        for (size_t t = m; t < m + rem; t++) {
            *q-- = bwt[ps];
            ps = lf[ps];
        }
    }
    free(lf);
    return 0;
}
