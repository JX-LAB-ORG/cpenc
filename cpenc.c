#include "cpenc.h"
#include "bwt.h"
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <stdio.h>

#define HUFF_MAX 12
#define HUFF_WSZ (1<<HUFF_MAX)
#define HUFF_TABSZ (HUFF_WSZ * 2)

/* ── entropy from freq table ─────────────────────────────────
 * fq_ent is called per group (tens of thousands per block); 256 libm
 * log2 calls each would dominate — precompute log2 for the integer
 * count range (segmented groups are <= 64KB so the table covers the
 * hot path; larger counts fall back to log2). */
static const uint8_t bitrev[256] = {
    0x00,0x80,0x40,0xC0,0x20,0xA0,0x60,0xE0,0x10,0x90,0x50,0xD0,0x30,0xB0,0x70,0xF0,
    0x08,0x88,0x48,0xC8,0x28,0xA8,0x68,0xE8,0x18,0x98,0x58,0xD8,0x38,0xB8,0x78,0xF8,
    0x04,0x84,0x44,0xC4,0x24,0xA4,0x64,0xE4,0x14,0x94,0x54,0xD4,0x34,0xB4,0x74,0xF4,
    0x0C,0x8C,0x4C,0xCC,0x2C,0xAC,0x6C,0xEC,0x1C,0x9C,0x5C,0xDC,0x3C,0xBC,0x7C,0xFC,
    0x02,0x82,0x42,0xC2,0x22,0xA2,0x62,0xE2,0x12,0x92,0x52,0xD2,0x32,0xB2,0x72,0xF2,
    0x0A,0x8A,0x4A,0xCA,0x2A,0xAA,0x6A,0xEA,0x1A,0x9A,0x5A,0xDA,0x3A,0xBA,0x7A,0xFA,
    0x06,0x86,0x46,0xC6,0x26,0xA6,0x66,0xE6,0x16,0x96,0x56,0xD6,0x36,0xB6,0x76,0xF6,
    0x0E,0x8E,0x4E,0xCE,0x2E,0xAE,0x6E,0xEE,0x1E,0x9E,0x5E,0xDE,0x3E,0xBE,0x7E,0xFE,
    0x01,0x81,0x41,0xC1,0x21,0xA1,0x61,0xE1,0x11,0x91,0x51,0xD1,0x31,0xB1,0x71,0xF1,
    0x09,0x89,0x49,0xC9,0x29,0xA9,0x69,0xE9,0x19,0x99,0x59,0xD9,0x39,0xB9,0x79,0xF9,
    0x05,0x85,0x45,0xC5,0x25,0xA5,0x65,0xE5,0x15,0x95,0x55,0xD5,0x35,0xB5,0x75,0xF5,
    0x0D,0x8D,0x4D,0xCD,0x2D,0xAD,0x6D,0xED,0x1D,0x9D,0x5D,0xDD,0x3D,0xBD,0x7D,0xFD,
    0x03,0x83,0x43,0xC3,0x23,0xA3,0x63,0xE3,0x13,0x93,0x53,0xD3,0x33,0xB3,0x73,0xF3,
    0x0B,0x8B,0x4B,0xCB,0x2B,0xAB,0x6B,0xEB,0x1B,0x9B,0x5B,0xDB,0x3B,0xBB,0x7B,0xFB,
    0x07,0x87,0x47,0xC7,0x27,0xA7,0x67,0xE7,0x17,0x97,0x57,0xD7,0x37,0xB7,0x77,0xF7,
    0x0F,0x8F,0x4F,0xCF,0x2F,0xAF,0x6F,0xEF,0x1F,0x9F,0x5F,0xDF,0x3F,0xBF,0x7F,0xFF};
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
    while(o--) for(size_t i=1;i<n;i++) p[i]+=p[i-1];
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
static void hlens(uint32_t fq[256],uint8_t lens[256]){
    memset(lens,0,256);
    HPair heap[512]; int n=0;
    for(int i=0;i<256;i++) if(fq[i]) hp_push(heap,&n,fq[i],i);
    if(n==0) return;
    if(n==1){lens[heap[0].sym]=1;return;}
    int left[512],right[512],next_id=256;
    memset(left,-1,sizeof(left)); memset(right,-1,sizeof(right));
    while(n>1){HPair a=hp_pop(heap,&n),b=hp_pop(heap,&n);left[next_id]=a.sym;right[next_id]=b.sym;hp_push(heap,&n,a.f+b.f,next_id++);}
    int root=heap[0].sym,q[512],qh=0,qt=0,depth[512];
    memset(depth,0,sizeof(depth));
    q[qt++]=root; depth[root]=0;
    /* Leaf depths are set here from their parent's depth; re-visiting a
     * leaf in the BFS would recompute d = depth[v]+1 and add 1 to every
     * symbol's code length (a long-standing -1 bit/symbol penalty). */
    while(qh<qt){int v=q[qh++];
        if(v>=256){int d=depth[v]+1;
            if(left[v]>=0){depth[left[v]]=d;q[qt++]=left[v];}
            if(right[v]>=0){depth[right[v]]=d;q[qt++]=right[v];}
        }
    }
    /* Length-limited adjustment (JPEG Annex K.3 style).  Clamping deep
     * leaves to HUFF_MAX would collide canonical codes; instead count
     * over-long leaves and rebalance shorter lengths to free code space. */
    int cnt[33]={0};
    for(int i=0;i<256;i++) if(depth[i]) cnt[depth[i]<=32?depth[i]:32]++;
    int excess=0;
    for(int d=HUFF_MAX+1;d<=32;d++) excess+=cnt[d];
    int li=HUFF_MAX-1;
    while(li>1&&excess>0){
        if(cnt[li]>1){cnt[li]-=2;cnt[li+1]++;excess--;}
        else li--;
    }
    /* Canonical re-assignment: most frequent symbols get the shortest codes */
    int order[256],no=0;
    for(int s=0;s<256;s++) if(fq[s]) order[no++]=s;
    for(int a=0;a<no;a++) for(int b=a+1;b<no;b++)
        if(fq[order[b]]>fq[order[a]]){int t=order[a];order[a]=order[b];order[b]=t;}
    int oi=0;
    for(int d=1;d<=HUFF_MAX&&oi<no;d++)
        for(int c=0;c<cnt[d]&&oi<no;c++) lens[order[oi++]]=(uint8_t)d;
    for(;oi<no;oi++) lens[order[oi]]=HUFF_MAX;
}
static void hcodes(const uint8_t lens[256],uint16_t code[256],uint8_t bits[256]){
    int cnt[16]={0};for(int i=0;i<256;i++)if(lens[i])cnt[lens[i]]++;
    uint32_t cur=0;uint16_t next[16];
    for(int L=1;L<=15;L++){next[L]=(uint16_t)cur;cur=(cur+(uint32_t)cnt[L])<<1;if(cur>65535)cur=65535;}
    for(int i=0;i<256;i++){if(lens[i]){code[i]=next[lens[i]]++;bits[i]=lens[i];}else{code[i]=0;bits[i]=0;}}
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
    if (getenv("CPENC_NO_THREADS")) return 1;
    if (!ts) return 1;
    int n = atoi(ts);
    if (n < 1) return 1;
    if (n > 16) n = 16;
    if ((size_t)n > n_seg) n = (int)n_seg;
    return n;
}
typedef struct { const uint8_t *src; size_t len; uint8_t *out; size_t cap; size_t outlen; int rc; int mo; } SegEncJob;
static void *seg_enc_worker(void *arg) {
    SegEncJob *j = arg;
    size_t dl = j->cap;
    j->rc = cpenc_compress(j->src, j->len, j->out, &dl, j->mo);
    j->outlen = dl;
    return NULL;
}
typedef struct { const uint8_t *src; size_t clen; uint8_t *out; size_t cap; size_t outlen; int rc; } SegDecJob;
static void *seg_dec_worker(void *arg) {
    SegDecJob *j = arg;
    size_t dl = j->cap;
    j->rc = cpenc_decompress(j->src, j->clen, j->out, &dl);
    j->outlen = dl;
    return NULL;
}
/* ── const+exceptions sub-block (subtag 4): [len u32][cval u8][nex u32]
 *    [pos LEB128][val u8] each — for streams dominated by one symbol ── */
static size_t sub5_estimate(const uint32_t cnt[256], size_t n) {
    uint32_t dom = 0;
    for (int i = 0; i < 256; i++) if (cnt[i] > dom) dom = cnt[i];
    return 9 + (size_t)(n - dom) * 5;
}
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
    if ((size_t)(edd - p) < 9) return -1;
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
#include <pthread.h>
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
    if ((size_t)(end - src) < 1 + 4*HUFF_SEGS) return -1;
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
    for (int k = 1; k < HUFF_SEGS; k++) segp[k] = segp[k-1] + lens[k-1];
    if ((size_t)(end - segp[HUFF_SEGS-1]) < lens[HUFF_SEGS-1]) return -1;
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
    const uint8_t *src=*sp;if((size_t)(end-src)<264)return-1;
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
    const uint8_t *src=*sp;if((size_t)(end-src)<264)return-1;
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

/* ═══════════════════════════════════════════════════════════════
   Public API
   ═══════════════════════════════════════════════════════════════ */

size_t cpenc_max_compressed_size(size_t sl) {
    return CPENC_HEADER_SIZE+16+sl+sl/2+16384;
}

int cpenc_compress_impl(const void *sv, size_t sl, void *dv, size_t *dl, int mo, int no_bwt) {
    lg2_init();
    const uint8_t *s=(const uint8_t*)sv; uint8_t *d=(uint8_t*)dv;
    CPEncHeader *h=(CPEncHeader*)d;
    if(*dl<CPENC_HEADER_SIZE) return -1;
    if(sl>0xFFFFFFFFU) return -2;
    mo=mo<0?0:mo>4?4:mo;
    if(sl==0){memset(h,0,CPENC_HEADER_SIZE);h->magic=CPENC_MAGIC;
        h->version=CPENC_VERSION;h->stride=CPENC_STRIDE_RAW;*dl=CPENC_HEADER_SIZE;return 0;}
    /* Tiny inputs: the stride/group machinery breaks down (gs can be 0)
     * and compression gains are nil — store RAW, always safe to decode. */
    if(sl<64){
        if(*dl<CPENC_HEADER_SIZE+1+sl) return -1;
        h->magic=CPENC_MAGIC;h->version=CPENC_VERSION;h->stride=CPENC_STRIDE_RAW;
        h->original_size=(uint32_t)sl;h->compressed_size=(uint32_t)(1+sl);
        d[CPENC_HEADER_SIZE]=CPENC_LAG_NONE;
        memcpy(d+CPENC_HEADER_SIZE+1,s,sl);*dl=CPENC_HEADER_SIZE+1+sl;return 0;
    }

    /* v4: hex-ascii layer — text-form numbers (log hex dumps) waste the
     * delta machinery on half-byte interleave ('F' 'F' = 0xFF, two chars
     * per value).  Detect hex-heavy text, fold 2 chars -> 1 byte, and keep
     * the non-hex layout chars (newlines etc.) in a position table so the
     * byte stream is restored exactly.  Encoded as [bin_len u32][n_sep u32]
     * [(idx LEB128)(char)]*n_sep [bin bytes]. */
    int hex_mode = 0;
    uint8_t *hexbuf = NULL;
    if (sl >= 256) {
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
            size_t blk = 13 + sl + sl / 2 + 16; /* hdr + table scratch + bin scratch */
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
                            hexbuf[13 + sl + n_sep++] = (uint8_t)(pend < 10 ? '0' + pend : 'a' + pend - 10);
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
                    hexbuf[13 + sl + n_sep++] = (uint8_t)(pend < 10 ? '0' + pend : 'a' + pend - 10);
                }
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
     * than the threshold recurse as plain single-block calls. */
    if (sl >= CPENC_SEG_MIN) {
        size_t n_seg = (sl + CPENC_SEG_SIZE - 1) / CPENC_SEG_SIZE;
        size_t ext = hex_mode ? 1 : 0;      /* real-stride byte for hex wrapper */
        size_t need = CPENC_HEADER_SIZE + ext + 4 + n_seg * 4;
        if (*dl < need) { free(hexbuf); return -1; }
        h->magic = CPENC_MAGIC; h->version = CPENC_VERSION;
        h->stride = hex_mode ? CPENC_STRIDE_HEX : CPENC_STRIDE_SEG;
        h->original_size = (uint32_t)sl;
        uint8_t *op = d + CPENC_HEADER_SIZE;
        if (hex_mode) *op++ = CPENC_STRIDE_SEG;
        *op++ = (uint8_t)(n_seg & 0xFF); *op++ = (uint8_t)((n_seg >> 8) & 0xFF);
        *op++ = (uint8_t)((n_seg >> 16) & 0xFF); *op++ = (uint8_t)((n_seg >> 24) & 0xFF);
        uint8_t *lenp = op;
        op += n_seg * 4;
        int nthr = seg_threads(n_seg);
        if (nthr > 1) {
            /* parallel: each segment compresses into its own buffer,
             * then a sequential pass assembles the stream */
            uint8_t **bufs = malloc(n_seg * sizeof(uint8_t*));
            SegEncJob *jobs = calloc(n_seg, sizeof(SegEncJob));
            if (!bufs || !jobs) { free(bufs); free(jobs); free(hexbuf); return -1; }
            for (size_t k = 0; k < n_seg; k++) {
                bufs[k] = malloc(cpenc_max_compressed_size(CPENC_SEG_SIZE) + 64);
                if (!bufs[k]) { for (size_t j = 0; j < k; j++) free(bufs[j]); free(bufs); free(jobs); free(hexbuf); return -1; }
                jobs[k].src = s + k * CPENC_SEG_SIZE;
                jobs[k].len = sl - k * CPENC_SEG_SIZE;
                if (jobs[k].len > CPENC_SEG_SIZE) jobs[k].len = CPENC_SEG_SIZE;
                jobs[k].out = bufs[k];
                jobs[k].cap = cpenc_max_compressed_size(CPENC_SEG_SIZE) + 64;
                jobs[k].mo = mo;
            }
            pthread_t th[16];
            for (size_t base = 0; base < n_seg; base += (size_t)nthr) {
                int nb = (int)(n_seg - base < (size_t)nthr ? n_seg - base : (size_t)nthr);
                for (int i = 0; i < nb; i++) pthread_create(&th[i], NULL, seg_enc_worker, &jobs[base + i]);
                for (int i = 0; i < nb; i++) pthread_join(th[i], NULL);
            }
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
                size_t seg_off = k * CPENC_SEG_SIZE;
                size_t seg_len = sl - seg_off;
                if (seg_len > CPENC_SEG_SIZE) seg_len = CPENC_SEG_SIZE;
                size_t clen = (size_t)(*dl - (op - d));
                if (cpenc_compress_impl(s + seg_off, seg_len, op, &clen, mo, 0) != 0) {
                    free(hexbuf); return -1;
                }
                lenp[k * 4 + 0] = (uint8_t)(clen & 0xFF);
                lenp[k * 4 + 1] = (uint8_t)((clen >> 8) & 0xFF);
                lenp[k * 4 + 2] = (uint8_t)((clen >> 16) & 0xFF);
                lenp[k * 4 + 3] = (uint8_t)((clen >> 24) & 0xFF);
                op += clen;
            }
        }
        h->compressed_size = (uint32_t)(op - d - CPENC_HEADER_SIZE);
        *dl = (size_t)(op - d);
        free(hexbuf);
        return 0;
    }

    /* v6: BWT candidate — reorder so context-similar bytes are adjacent,
     * then the lag/stride machinery eats the order-1 structure ("BWT is
     * the delta aligner").  Archive-oriented: sort is slow (~3MB/s), the
     * LF-mapping decode is O(n).  Mid-entropy blocks only; structured
     * (<4.5) and random (>7.2) data skip — entropy gate is the cost cap.
     * Format: [orig_idx u32][payload] where payload = compressed BWT stream
     * (recursed with no_bwt=1 so it never re-enters this candidate). */
    if (!no_bwt && !hex_mode && sl >= 8192 && sl <= (1u<<20)) {
        double e = ent(s, 2048);
        if (e > 4.5 && e < 7.2) {
            uint8_t *bw = malloc(sl + 1);
            if (bw) {
                size_t oi = bwt_enc(s, sl, bw);
                /* bzip2-style: MTF on the BWT stream turns its local
                 * clustering into long zero runs (RLE territory) */
                mtf_enc_inplace(bw, sl + 1);
                uint32_t bc[256] = {0};
                for (size_t i = 0; i <= sl; i++) bc[bw[i]]++;
                double e1 = fq_ent(bc, sl + 1);
                if (e1 < e - 0.3 && *dl > CPENC_HEADER_SIZE + 4 + 64) {
                    size_t cl = *dl - (CPENC_HEADER_SIZE + 4);
                    uint8_t *op = d + CPENC_HEADER_SIZE + 4;
                    if (cpenc_compress_impl(bw, sl + 1, op, &cl, mo, 1) == 0 &&
                        cl + CPENC_HEADER_SIZE + 4 < (size_t)(e * sl / 8.0 * 1.05 + 96)) {
                        h->magic = CPENC_MAGIC; h->version = CPENC_VERSION;
                        h->stride = CPENC_STRIDE_BWT;
                        h->original_size = (uint32_t)sl;
                        h->compressed_size = (uint32_t)(cl + 4);
                        d[CPENC_HEADER_SIZE+0] = (uint8_t)oi;
                        d[CPENC_HEADER_SIZE+1] = (uint8_t)(oi>>8);
                        d[CPENC_HEADER_SIZE+2] = (uint8_t)(oi>>16);
                        d[CPENC_HEADER_SIZE+3] = (uint8_t)(oi>>24);
                        *dl = CPENC_HEADER_SIZE + 4 + cl;
                        free(bw);
                        return 0;
                    }
                }
                free(bw);
            }
        }
    }

    /* v3: lag-delta preprocessing — catches fixed-frame periodicity
     * (log rows, protocol frames, config tables).  Coarse scan (every 16)
     * then fine scan around the best; high-entropy segments skip entirely
     * (no frame periodicity to find).  Diff + entropy merged into one pass. */
    int lag_idx = CPENC_LAG_NONE;
    uint8_t *lagbuf = NULL;
    if (sl >= 4096) {
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
        if (seg_e > 7.5) {
            size_t ext = hex_mode ? 2 : 1;
            if (*dl < CPENC_HEADER_SIZE + ext + sl) { free(lagbuf); free(hexbuf); return -1; }
            h->magic=CPENC_MAGIC;h->version=CPENC_VERSION;
            h->stride= hex_mode ? CPENC_STRIDE_HEX : CPENC_STRIDE_RAW;
            h->original_size=(uint32_t)sl;h->compressed_size=(uint32_t)(ext+sl);
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
    uint8_t *trial = malloc(4096); if(!trial){free(trans);return-2;}

    if(stride==1){memcpy(trans,s,sl);}
    else{tr(s,trans,sl,stride);}

    uint8_t orders[8];
    int nraw=0;

    for(int g=0;g<stride;g++){
        uint8_t *gp = trans + (size_t)g*gs;
        int best_o=0;
        size_t sample = gs<4096?gs:4096;
        double sample_e = ent(gp, sample);
        memcpy(trial,gp,sample); de(trial,sample,1);
        double e1=ent(trial,sample);
        if(e1<sample_e-0.3){
            best_o=1; double best_e=e1;
            if(mo>=2 && best_e>1.0){memcpy(trial,gp,sample);de(trial,sample,2);
                double e2=ent(trial,sample);if(e2<best_e-0.01){best_o=2;best_e=e2;}}
            if(mo>=3 && best_e>1.0){memcpy(trial,gp,sample);de(trial,sample,3);
                double e3=ent(trial,sample);if(e3<best_e-0.01){best_o=3;}}
        }else if(sample_e>6.5){best_o=0;nraw++;}
        de(gp,gs,best_o);
        orders[g]=(uint8_t)best_o;
    }

    /* Global RAW only for stride>1 where all groups independently
     * failed delta — a single-group stride=1 with high entropy may
     * still benefit from Huffman/LUT. Let per-group tag selection decide. */
    if(nraw==stride && stride>1 && sl>=256){
        size_t ext = hex_mode ? 2 : 1;
        if(*dl<CPENC_HEADER_SIZE+ext+sl){free(trial);free(trans);free(lagbuf);free(hexbuf);return-1;}
        h->magic=CPENC_MAGIC;h->version=CPENC_VERSION;h->stride= hex_mode ? CPENC_STRIDE_HEX : CPENC_STRIDE_RAW;
        h->original_size=(uint32_t)sl;h->compressed_size=(uint32_t)(ext+sl);
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
        double gp_ent=fq_ent(fq,gs);
        /* huffman min length is 1 bit/symbol for the dominant symbol:
         * the bitstream floor is gs/8 bytes, not gs/32 (was badly
         * underestimating skewed groups, steering them to Huffman) */
        size_t h_p=(size_t)(gs*gp_ent/8.0);if(h_p<gs/8)h_p=gs/8;
        /* zero-dense groups (BWT+MTF output, sparse): Huffman's actual code
         * length runs ~1.75x the entropy floor — the 0 gets 1 bit but the
         * long tail of sparse ranks eats the rest.  Without this, Huffman
         * estimate beats RLE and loses ~2x ratio (measured on MTF streams). */
        if (fq[0] > (uint32_t)gs/2) h_p = h_p * 175 / 100;
        size_t sz_h=264+h_p*103/100+16;

        if(nz<=1){tag=0;}
        else{
            int bits=1;while((1<<bits)<nz&&bits<HUFF_MAX)bits++;
            size_t sz_l=6+(size_t)nz+gs*(size_t)bits/8+5;
            size_t sz_best=sz_raw; tag=3;
            if(sz_l<sz_best){tag=1;sz_best=sz_l;}
            if(sz_h<sz_best){tag=2;sz_best=sz_h;}
            /* tag 5: single-value dominance — one symbol covers >99% of the
             * group.  Encode the constant plus a short exception list
             * (positions LEB128 + values), vs Huffman's wasteful 1-bit-per-
             * symbol all-zero bitstream.  Complements RLE (zeros) with
             * non-zero constants (e.g. transposed struct groups). */
            {
                uint32_t dom = 0;
                for (int i2 = 0; i2 < 256; i2++) if (fq[i2] > dom) { dom = fq[i2]; }
                if (dom > (uint32_t)gs - gs/10 && dom < (uint32_t)gs && gs >= 1024) {
                    /* estimate: cval(2) + count(4) + gap-LEB128 pos (avg 2B) + value */
                    uint32_t nex = (uint32_t)gs - dom;
                    size_t sz_c = 6 + (size_t)nex*3;
                    if (sz_c < sz_best) { tag = 5; }
                }
            }
            /* RLE candidate (tag 4): zero-heavy groups get run-length coded.
             * Strict alternation: zero-run then value; consecutive nonzeros
             * are separated by run=0.  Runs are 8-bit with 255 = continuation.
             * Threshold lowered to 50% so BWT+MTF streams (0 ~57-70%, long
             * runs) qualify; run count is scanned (cheap, only for zero-heavy
             * groups) because the alternating-stream bound 2nv+1 overestimates
             * clustered zeros by ~5x and would steer them to Huffman. */
            if(fq[0] > (uint32_t)gs/2 && gs >= 64){
                size_t nv = (size_t)gs - fq[0];
                size_t nr = 0; int inz = 0;
                for (size_t i2 = 0; i2 < gs; i2++) {
                    if (!gp[i2]) { if (!inz) { nr++; inz = 1; } }
                    else inz = 0;
                }
                /* value sub-stream gets its own tag; entropy of the nonzero
                 * symbols (RAW sizing overestimated it ~40% on MTF streams,
                 * steering RLE to Huffman) */
                double ve = 0;
                if (nv) {
                    for (int i2 = 1; i2 < 256; i2++)
                        if (fq[i2]) { double pr = (double)fq[i2]/nv; ve -= pr*log2(pr); }
                }
                size_t v_est = 264 + (size_t)(nv*ve/8.0*103/100) + 16;
                size_t sz_rle = 2 + (4 + nr + nr/5) + (4 + v_est);
                if(sz_rle<sz_best){tag=4;}
            }
        }

        if(out>=edd){free(trial);free(trans);free(lagbuf);free(hexbuf);return-1;} *out++=tag;
        if(tag==4){
            /* RLE: rebuild streams and encode both sub-blocks with the
             * standard per-subtag machinery (constants/LUT/Huffman/RAW) */
            uint32_t *runs = malloc((gs/2+2)*sizeof(uint32_t));
            uint8_t *vls = malloc(gs);
            if(!runs||!vls){free(runs);free(vls);free(trial);free(trans);free(lagbuf);free(hexbuf);return-1;}
            size_t nr=0,nv=0,i;
            if(gp[0]!=0) runs[nr++]=0;
            for(i=0;i<gs;){
                if(gp[i]==0){size_t j=i;while(j<gs&&gp[j]==0)j++;
                    runs[nr++]=(uint32_t)(j-i);i=j;}
                else{if(i>0&&gp[i-1]!=0)runs[nr++]=0;
                    vls[nv++]=gp[i];i++;}
            }
            uint8_t *r8=malloc(nr*5+2);
            if(!r8){free(runs);free(vls);free(trial);free(trans);free(lagbuf);free(hexbuf);return-1;}
            size_t nr8=0;
            for(size_t k=0;k<nr;k++){
                if(runs[k]<255){r8[nr8++]=(uint8_t)runs[k];}
                else{r8[nr8++]=255;   /* escape: length-255 remainder in LEB128 */
                    uint32_t r=runs[k]-255;
                    do{uint8_t b=(uint8_t)(r&0x7F);r>>=7;
                        if(r)b|=0x80;
                        r8[nr8++]=b;}while(r);}
            }

            /* encode runs sub-block (subtag) */
            uint32_t rq[256]={0};for(size_t k=0;k<nr8;k++)rq[r8[k]]++;
            int rnz=0;for(int i2=0;i2<256;i2++)if(rq[i2])rnz++;
            uint8_t st;
            if(rnz<=1)st=0;
            else{
                int rb=1;while((1<<rb)<rnz&&rb<HUFF_MAX)rb++;
                size_t sl=6+rnz+nr8*(size_t)rb/8+5;
                double re=fq_ent(rq,nr8);size_t hp=(size_t)(nr8*re/8.0);if(hp<nr8/8)hp=nr8/8;
                size_t sh=264+hp*103/100+16;
                size_t raw=4+nr8;
                size_t s5=sub5_estimate(rq,nr8);
                size_t best=sl; st=1;
                if(sh<best){best=sh;st=2;}
                if(raw<best){best=raw;st=3;}
                if(s5<best){best=s5;st=4;}
            }
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
                uint8_t lut[32],rev[256]={0};int li=0;
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
            if(vnz<=1)st=0;
            else{
                int vb=1;while((1<<vb)<vnz&&vb<HUFF_MAX)vb++;
                size_t sl=6+vnz+nv*(size_t)vb/8+5;
                double ve=fq_ent(vq,nv);size_t hp=(size_t)(nv*ve/8.0);if(hp<nv/8)hp=nv/8;
                size_t sh=264+hp*103/100+16;
                size_t raw=4+nv;
                size_t s5=sub5_estimate(vq,nv);
                size_t best=sl; st=1;
                if(sh<best){best=sh;st=2;}
                if(raw<best){best=raw;st=3;}
                if(s5<best){best=s5;st=4;}
            }
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
                uint8_t lut[32],rev[256]={0};int li=0;
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
            uint8_t lut[32],rev[256]={0};int li=0;
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
        }else{
            if((size_t)(edd-out)<4+gs){free(trial);free(trans);free(lagbuf);free(hexbuf);return-1;}
            *out++=gs&0xFF;*out++=(gs>>8)&0xFF;*out++=(gs>>16)&0xFF;*out++=gs>>24;
            memcpy(out,gp,gs);out+=gs;
        }
    }

    h->magic=CPENC_MAGIC;h->version=CPENC_VERSION;h->stride= hex_mode ? CPENC_STRIDE_HEX : (uint8_t)stride;
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
    h->original_size=(uint32_t)sl;h->compressed_size=(uint32_t)(out-(d+CPENC_HEADER_SIZE));
    *dl=(size_t)(out-d); free(trial); free(trans); free(hexbuf); return 0;
}

int cpenc_decompress(const void *sv, size_t sl, void *dv, size_t *dl) {
    const uint8_t *s=(const uint8_t*)sv; uint8_t *d=(uint8_t*)dv;
    const CPEncHeader *h=(const CPEncHeader*)s;
    if(sl<CPENC_HEADER_SIZE) return -1;
    /* v1 blocks (no stride remainder) remain decodable; reject only
     * unknown versions and garbage magic */
    if(h->magic!=CPENC_MAGIC||h->version<1||h->version>CPENC_VERSION) return -1;
    if(h->stride==CPENC_STRIDE_SEG){
        /* v5: independent sub-blocks, lengths in a u32 array */
        const uint8_t *p = s + CPENC_HEADER_SIZE;
        if(sl < CPENC_HEADER_SIZE + 4) return -1;
        uint32_t n_seg = (uint32_t)p[0] | ((uint32_t)p[1]<<8) | ((uint32_t)p[2]<<16) | ((uint32_t)p[3]<<24);
        p += 4;
        uint32_t max_seg = (h->original_size + CPENC_SEG_SIZE - 1) / CPENC_SEG_SIZE;
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
            if (!bufs || !jobs) { free(bufs); free(jobs); return -1; }
            for (uint32_t k = 0; k < n_seg; k++) {
                uint32_t clen = (uint32_t)p[k*4] | ((uint32_t)p[k*4+1]<<8) | ((uint32_t)p[k*4+2]<<16) | ((uint32_t)p[k*4+3]<<24);
                size_t dl_k = h->original_size - (size_t)k * CPENC_SEG_SIZE;
                if (dl_k > CPENC_SEG_SIZE) dl_k = CPENC_SEG_SIZE;
                bufs[k] = malloc(dl_k ? dl_k : 1);
                if (!bufs[k]) { for (uint32_t j = 0; j < k; j++) free(bufs[j]); free(bufs); free(jobs); return -1; }
                jobs[k].src = sp; jobs[k].clen = clen;
                jobs[k].out = bufs[k]; jobs[k].cap = dl_k;
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
                if (jobs[k].rc != 0 || jobs[k].outlen != jobs[k].cap) {
                    for (uint32_t j = 0; j < n_seg; j++) free(bufs[j]);
                    free(bufs); free(jobs); return -1;
                }
                if ((size_t)*dl < o + jobs[k].outlen) {
                    for (uint32_t j = 0; j < n_seg; j++) free(bufs[j]);
                    free(bufs); free(jobs); return -1;
                }
                memcpy(d + o, bufs[k], jobs[k].outlen);
                o += jobs[k].outlen;
            }
            for (uint32_t j = 0; j < n_seg; j++) free(bufs[j]);
            free(bufs); free(jobs);
            *dl = o;
            return 0;
        }
        size_t o = 0;
        for(uint32_t k = 0; k < n_seg; k++) {
            uint32_t clen = (uint32_t)p[k*4] | ((uint32_t)p[k*4+1]<<8) | ((uint32_t)p[k*4+2]<<16) | ((uint32_t)p[k*4+3]<<24);
            size_t dl_k = h->original_size - o;
            if(dl_k > CPENC_SEG_SIZE) dl_k = CPENC_SEG_SIZE;
            if((size_t)*dl < o + dl_k) return -1;
            if(cpenc_decompress(sp, clen, d + o, &dl_k) != 0) return -1;
            size_t want = h->original_size - o;
            if(want > CPENC_SEG_SIZE) want = CPENC_SEG_SIZE;
            if(dl_k != want) return -1;
            o += dl_k;
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
                for(;;){ uint8_t b = tab[oi++]; idx |= (uint32_t)(b & 0x7F) << sh;
                         if(!(b & 0x80)) break; sh += 7;
                         if(oi >= n_sep && (b & 0x80)){ free(tmp); return -1; } }
                if(idx == vi) { d[o++] = tab[oi++];  }
                else { oi = save; break; }
            }
            d[o++] = (uint8_t)hexd[bin[i] >> 4];
            d[o++] = (uint8_t)hexd[bin[i] & 0xF];
            vi++;
        }
        while(oi < n_sep) {                  /* trailing entries (idx == bin_len) */
            uint32_t idx = 0; int sh = 0;
            for(;;){ uint8_t b = tab[oi++]; idx |= (uint32_t)(b & 0x7F) << sh;
                     if(!(b & 0x80)) break; sh += 7; }
            if(idx != vi){ free(tmp); return -1; }
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
        uint8_t *bw = malloc(h->original_size + 1);
        if(!bw) return -2;
        size_t outl = h->original_size + 1;
        int rc = cpenc_decompress(rp + 4, pl, bw, &outl);
        if(rc != 0 || outl != h->original_size + 1){ free(bw); return -1; }
        mtf_dec_inplace(bw, outl);               /* undo MTF -> BWT stream */
        uint8_t *tmp = malloc(outl);
        if(!tmp){ free(bw); return -2; }
        bwt_dec(bw, outl, oi, tmp);              /* n+1 rows, EOF 0 first */
        memcpy(d, tmp + 1, h->original_size);    /* strip the EOF symbol */
        free(tmp);
        free(bw);
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

    for(int g=0;g<stride;g++){
        if(p>=edd){if(stride>1)free(trans);return-1;}
        uint8_t tag=*p++;
        uint8_t *gp=trans+g*gs;

        if(tag==3){
            if((size_t)(edd-p)<4){if(stride>1)free(trans);return-1;}
            uint32_t rl=(uint32_t)p[0]|((uint32_t)p[1]<<8)|((uint32_t)p[2]<<16)|((uint32_t)p[3]<<24);p+=4;
            if(rl!=gs||(size_t)(edd-p)<rl){if(stride>1)free(trans);return-1;}
            memcpy(gp,p,rl);p+=rl;
        }else if(tag==0){
            if((size_t)(edd-p)<6){if(stride>1)free(trans);return-1;}
            uint32_t tgs=(uint32_t)p[0]|((uint32_t)p[1]<<8)|((uint32_t)p[2]<<16)|((uint32_t)p[3]<<24);p+=4;
            uint8_t val=*p++;p++;
            if(tgs!=gs){if(stride>1)free(trans);return-1;}
            memset(gp,val,gs);
        }else if(tag==1){
            if((size_t)(edd-p)<6){if(stride>1)free(trans);return-1;}
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
        }else if(tag==5){
            /* const + exceptions: [val][nex][pos LEB128][values] */
            if((size_t)(edd-p)<6){if(stride>1)free(trans);return-1;}
            uint16_t cval=(uint16_t)(p[0]|(p[1]<<8));p+=2;
            uint32_t nex=(uint32_t)p[0]|((uint32_t)p[1]<<8)|((uint32_t)p[2]<<16)|((uint32_t)p[3]<<24);p+=4;
            if(nex>gs){if(stride>1)free(trans);return-1;}
            memset(gp,cval,gs);
            uint32_t pos=0;
            for(uint32_t k=0;k<nex;k++){
                uint32_t gap=0,shift=0;
                for(;;){
                    if(p>=edd){if(stride>1)free(trans);return-1;}
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
            uint32_t *runs_buf = malloc((gs/2+2)*sizeof(uint32_t));
            uint8_t *vals_buf = malloc(gs);
            uint8_t *scratch = malloc(gs);
            if(!runs_buf||!vals_buf||!scratch){free(runs_buf);free(vals_buf);free(scratch);if(stride>1)free(trans);return-1;}
            size_t nr=0,nv=0;
            if(st==0){
                if((size_t)(edd-p)<6){free(runs_buf);free(vals_buf);free(scratch);if(stride>1)free(trans);return-1;}
                uint32_t tgs=(uint32_t)p[0]|((uint32_t)p[1]<<8)|((uint32_t)p[2]<<16)|((uint32_t)p[3]<<24);p+=4;
                uint8_t val=*p++;p++;
                if(tgs>gs/2+1){free(runs_buf);free(vals_buf);free(scratch);if(stride>1)free(trans);return-1;}
                for(size_t k=0;k<tgs;k++)runs_buf[nr++]=(uint32_t)val;
            }else if(st==1){
                if((size_t)(edd-p)<6){free(runs_buf);free(vals_buf);free(scratch);if(stride>1)free(trans);return-1;}
                uint32_t tgs=(uint32_t)p[0]|((uint32_t)p[1]<<8)|((uint32_t)p[2]<<16)|((uint32_t)p[3]<<24);p+=4;
                uint8_t nz1=*p++;int nz=nz1+1;
                uint8_t lut[256];for(int i=0;i<nz;i++){if(p>=edd){free(runs_buf);free(vals_buf);free(scratch);if(stride>1)free(trans);return-1;}lut[i]=*p++;}
                if(tgs>gs/2+1){free(runs_buf);free(vals_buf);free(scratch);if(stride>1)free(trans);return-1;}
                int nbits=1;while((1<<nbits)<nz&&nbits<12)nbits++;
                uint64_t rbuf=0;int rbc=0;
                for(size_t i=0;i<tgs;i++){
                    while(rbc<nbits&&p<edd){rbuf|=(uint64_t)(*p++)<<rbc;rbc+=8;}
                    uint8_t idx=rbuf&((1<<nbits)-1);rbuf>>=nbits;rbc-=nbits;
                    runs_buf[nr++]=(uint32_t)(idx<nz?lut[idx]:lut[0]);
                }
            }else if(st==2){
                /* huffman sub-block: symbol count lives in the stream header */
                if((size_t)(edd-p)<264){free(runs_buf);free(vals_buf);free(scratch);if(stride>1)free(trans);return-1;}
                uint32_t rn=(uint32_t)p[256]|((uint32_t)p[257]<<8)|((uint32_t)p[258]<<16)|((uint32_t)p[259]<<24);
                if(rn>gs/2+2){free(runs_buf);free(vals_buf);free(scratch);if(stride>1)free(trans);return-1;}
                if(h_dec(&p,edd,scratch,rn)<0){free(runs_buf);free(vals_buf);free(scratch);if(stride>1)free(trans);return-1;}
                for(size_t k=0;k<rn;k++)runs_buf[nr++]=(uint32_t)scratch[k];
            }else if(st==3){
                if((size_t)(edd-p)<4){free(runs_buf);free(vals_buf);free(scratch);if(stride>1)free(trans);return-1;}
                uint32_t rl=(uint32_t)p[0]|((uint32_t)p[1]<<8)|((uint32_t)p[2]<<16)|((uint32_t)p[3]<<24);p+=4;
                if(rl>gs/2+1||(size_t)(edd-p)<rl){free(runs_buf);free(vals_buf);free(scratch);if(stride>1)free(trans);return-1;}
                for(size_t k=0;k<rl;k++)runs_buf[nr++]=(uint32_t)p[k];
                p+=rl;
            }else if(st==4){
                uint32_t tlen=0;
                if(sub5_decode(&p,edd,scratch,gs/2+2,&tlen)!=0){free(runs_buf);free(vals_buf);free(scratch);if(stride>1)free(trans);return-1;}
                for(size_t k=0;k<tlen;k++)runs_buf[nr++]=(uint32_t)scratch[k];
            }
            /* vals sub-block */

            if(p>=edd){free(runs_buf);free(vals_buf);free(scratch);if(stride>1)free(trans);return-1;}
            st=*p++;

            if(st==0){
                if((size_t)(edd-p)<6){free(runs_buf);free(vals_buf);free(scratch);if(stride>1)free(trans);return-1;}
                uint32_t tgs=(uint32_t)p[0]|((uint32_t)p[1]<<8)|((uint32_t)p[2]<<16)|((uint32_t)p[3]<<24);p+=4;
                uint8_t val=*p++;p++;
                if(tgs>gs){free(runs_buf);free(vals_buf);free(scratch);if(stride>1)free(trans);return-1;}
                for(size_t k=0;k<tgs;k++)vals_buf[nv++]=(uint8_t)val;
            }else if(st==1){
                if((size_t)(edd-p)<6){free(runs_buf);free(vals_buf);free(scratch);if(stride>1)free(trans);return-1;}
                uint32_t tgs=(uint32_t)p[0]|((uint32_t)p[1]<<8)|((uint32_t)p[2]<<16)|((uint32_t)p[3]<<24);p+=4;
                uint8_t nz1=*p++;int nz=nz1+1;
                uint8_t lut[256];for(int i=0;i<nz;i++){if(p>=edd){free(runs_buf);free(vals_buf);free(scratch);if(stride>1)free(trans);return-1;}lut[i]=*p++;}
                if(tgs>gs){free(runs_buf);free(vals_buf);free(scratch);if(stride>1)free(trans);return-1;}
                int nbits=1;while((1<<nbits)<nz&&nbits<12)nbits++;
                uint64_t rbuf=0;int rbc=0;
                for(size_t i=0;i<tgs;i++){
                    while(rbc<nbits&&p<edd){rbuf|=(uint64_t)(*p++)<<rbc;rbc+=8;}
                    uint8_t idx=rbuf&((1<<nbits)-1);rbuf>>=nbits;rbc-=nbits;
                    vals_buf[nv++]=(uint8_t)(idx<nz?lut[idx]:lut[0]);
                }
            }else if(st==2){
                if((size_t)(edd-p)<264){free(runs_buf);free(vals_buf);free(scratch);if(stride>1)free(trans);return-1;}
                uint32_t rn=(uint32_t)p[256]|((uint32_t)p[257]<<8)|((uint32_t)p[258]<<16)|((uint32_t)p[259]<<24);
                if(rn>gs){free(runs_buf);free(vals_buf);free(scratch);if(stride>1)free(trans);return-1;}
                if(h_dec(&p,edd,vals_buf,rn)<0){free(runs_buf);free(vals_buf);free(scratch);if(stride>1)free(trans);return-1;}
                nv=rn;
            }else if(st==3){
                if((size_t)(edd-p)<4){free(runs_buf);free(vals_buf);free(scratch);if(stride>1)free(trans);return-1;}
                uint32_t rl=(uint32_t)p[0]|((uint32_t)p[1]<<8)|((uint32_t)p[2]<<16)|((uint32_t)p[3]<<24);p+=4;
                if(rl>gs||(size_t)(edd-p)<rl){free(runs_buf);free(vals_buf);free(scratch);if(stride>1)free(trans);return-1;}
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
                        if(ri>=nr){free(runs_buf);free(vals_buf);free(scratch);if(stride>1)free(trans);return-1;}
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
            if(stride==1 && orders[g]==1 && lag_idx==0)
                hrc = h_dec_p(&p,edd,gp,gs,1);       /* fused prefix sum */
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
            if ((size_t)(edd - p) < rem) { if (stride > 1) free(trans); return -1; }
            memcpy(trans + gend, p, rem);
            p += rem;
        }
    }
    for(int g=0;g<stride;g++) dd(trans+g*gs,gs,orders[g]);
    if(stride>1){ur(trans,d,orig,stride);free(trans);}
    /* v3: undo lag-delta (applied pre-transpose on the raw byte order) */
    if(lag_idx){
        for(size_t i=(size_t)lag_idx;i<orig;i++) d[i]+=d[i-(size_t)lag_idx];
    }
    *dl=orig; return 0;
}

int cpenc_compress(const void *sv, size_t sl, void *dv, size_t *dl, int mo) {
    return cpenc_compress_impl(sv, sl, dv, dl, mo, 0);
}
