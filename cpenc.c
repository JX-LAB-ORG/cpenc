#include "cpenc.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#define HUFF_MAX 12
#define HUFF_WSZ (1<<HUFF_MAX)
#define HUFF_TABSZ (HUFF_WSZ * 2)

/* ── entropy from freq table ───────────────────────────────── */
static double fq_ent(const uint32_t fq[256], size_t tot) {
    if(!tot) return 0;
    double h=0, iv=1.0/(double)tot;
    for(int i=0;i<256;i++) if(fq[i]){double pr=(double)fq[i]*iv; h-=pr*log2(pr);}
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
    while(qh<qt){int v=q[qh++],d=depth[v]+1;
        if(v>=256){if(left[v]>=0){depth[left[v]]=d;q[qt++]=left[v];}if(right[v]>=0){depth[right[v]]=d;q[qt++]=right[v];}}
        else{depth[v]=d;}
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
    for(size_t i=0;i<n;i++){int s=d[i];buf=(buf<<bits[s])|code[s];cnt+=bits[s];
        while(cnt>=8){cnt-=8;*op++=(uint8_t)(buf>>cnt);buf&=((1ULL<<cnt)-1);}}
    if(cnt>0)*op++=(uint8_t)(buf<<(8-cnt));
    uint32_t bsl=(uint32_t)(op-bst);ob[260]=bsl&0xFF;ob[261]=(bsl>>8)&0xFF;ob[262]=(bsl>>16)&0xFF;ob[263]=bsl>>24;
    return (uint32_t)(op-ob);
}
static int h_dec(const uint8_t **sp,const uint8_t *end,uint8_t *dst,size_t n){
    const uint8_t *src=*sp;if((size_t)(end-src)<264)return-1;
    uint8_t lens[256];for(int i=0;i<256;i++)lens[i]=src[i];
    uint32_t ul=(uint32_t)src[256]|((uint32_t)src[257]<<8)|((uint32_t)src[258]<<16)|((uint32_t)src[259]<<24);
    uint32_t bsl=(uint32_t)src[260]|((uint32_t)src[261]<<8)|((uint32_t)src[262]<<16)|((uint32_t)src[263]<<24);
    src+=264;const uint8_t *bse=src+bsl;
    uint16_t code[256];uint8_t bits[256];hcodes(lens,code,bits);
    uint8_t fast[HUFF_TABSZ];memset(fast,0,sizeof(fast));
    for(int s=0;s<256;s++){if(lens[s]&&lens[s]<=HUFF_MAX){uint16_t c=code[s];int b=lens[s];
        int xm=1<<(HUFF_MAX-b);
        for(int x=0;x<xm;x++){uint32_t fi=((uint32_t)c<<(HUFF_MAX-b))|(uint32_t)x;if(fi<HUFF_TABSZ)fast[fi]=(uint8_t)s;}}}
    size_t tgt=ul<n?ul:n;uint64_t buf=0;int cnt=0;
    for(size_t dc=0;dc<tgt;dc++){while(cnt<HUFF_MAX&&src<bse){buf=(buf<<8)|(*src++);cnt+=8;}
        if(cnt<HUFF_MAX){buf<<=(HUFF_MAX-cnt);cnt=HUFF_MAX;}
        uint32_t idx=(uint32_t)((buf>>(cnt-HUFF_MAX))&(HUFF_WSZ-1));uint8_t s=fast[idx];
        if(lens[s]==0) return -1;
        dst[dc]=s; cnt-=lens[s];
    }
    *sp=bse;return 0;
}

/* ═══════════════════════════════════════════════════════════════
   Public API
   ═══════════════════════════════════════════════════════════════ */

size_t cpenc_max_compressed_size(size_t sl) {
    return CPENC_HEADER_SIZE+16+sl+sl/2+16384;
}

int cpenc_compress(const void *sv, size_t sl, void *dv, size_t *dl, int mo) {
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
        if(*dl<CPENC_HEADER_SIZE+sl) return -1;
        h->magic=CPENC_MAGIC;h->version=CPENC_VERSION;h->stride=CPENC_STRIDE_RAW;
        h->original_size=(uint32_t)sl;h->compressed_size=(uint32_t)sl;
        memcpy(d+CPENC_HEADER_SIZE,s,sl);*dl=CPENC_HEADER_SIZE+sl;return 0;
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
        size_t gs_test=sl/(size_t)stride;
        uint8_t *tt=malloc(sl),*tb=malloc(gs_test);
        int nraw=0;
        if(tt&&tb){
            tr(s,tt,sl,stride);
            for(int g=0;g<stride;g++){uint8_t *gp=tt+(size_t)g*gs_test;
                memcpy(tb,gp,gs_test);double e0=ent(tb,gs_test);
                de(tb,gs_test,1);double e1=ent(tb,gs_test);
                if(e0>6.5&&e1>6.5)nraw++;}
            if(nraw==stride)stride=1;
        }
        free(tt);free(tb);
    }

    size_t gs=sl/(size_t)stride;
    uint8_t *trans = malloc(sl); if(!trans) return -2;
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
        if(*dl<CPENC_HEADER_SIZE+sl){free(trial);free(trans);return-1;}
        h->magic=CPENC_MAGIC;h->version=CPENC_VERSION;h->stride=CPENC_STRIDE_RAW;
        h->original_size=(uint32_t)sl;h->compressed_size=(uint32_t)sl;
        memcpy(d+CPENC_HEADER_SIZE,s,sl);*dl=CPENC_HEADER_SIZE+sl;
        free(trial);free(trans);return 0;
    }

    uint8_t *out=d+CPENC_HEADER_SIZE,*edd=d+*dl;
    for(int g=0;g<stride;g++){if(out>=edd){free(trial);free(trans);return-1;}*out++=orders[g];}

    for(int g=0;g<stride;g++){
        uint8_t *gp = trans + (size_t)g*gs;
        uint32_t fq[256];
        freq256(gp, gs, fq);
        int nz=0;for(int i=0;i<256;i++)if(fq[i])nz++;

        size_t sz_raw=4+gs; uint8_t tag;
        double gp_ent=fq_ent(fq,gs);
        size_t h_p=(size_t)(gs*gp_ent/8.0);if(h_p<gs/32)h_p=gs/32;
        size_t sz_h=264+h_p*103/100+16;

        if(nz<=1){tag=0;}
        else if(nz<=32){
            int bits=1;while((1<<bits)<nz&&bits<HUFF_MAX)bits++;
            size_t sz_l=6+(size_t)nz+gs*(size_t)bits/8+5;
            tag = sz_l<sz_raw ? 1 : 3;
        }
        else{tag=sz_h<sz_raw?2:3;}

        if(out>=edd){free(trial);free(trans);return-1;} *out++=tag;
        if(tag==0){
            uint8_t val=0;for(int i=0;i<256;i++)if(fq[i]){val=(uint8_t)i;break;}
            if(out+6>edd){free(trial);free(trans);return-1;}
            *out++=gs&0xFF;*out++=(gs>>8)&0xFF;*out++=(gs>>16)&0xFF;*out++=gs>>24;
            *out++=val;*out++=(uint8_t)(nz-1);
        }else if(tag==1){
            int bits=1;while((1<<bits)<nz&&bits<HUFF_MAX)bits++;
            uint8_t lut[32],rev[256]={0};int li=0;
            for(int i=0;i<256&&li<nz;i++)if(fq[i]){uint8_t v=(uint8_t)i;lut[li]=v;rev[v]=(uint8_t)li;li++;}
            size_t need=6+(size_t)nz+gs*(size_t)bits/8+1;
            if((size_t)(edd-out)<need){free(trial);free(trans);return-1;}
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
            uint8_t lens[256];uint16_t code[256];uint8_t bits[256];
            hlens(fq,lens);hcodes(lens,code,bits);
            uint32_t used=h_enc(gp,gs,lens,code,bits,out,(size_t)(edd-out));
            if(used==0||(size_t)(edd-out)<used){free(trial);free(trans);return-1;}
            out+=used;
        }else{
            if((size_t)(edd-out)<4+gs){free(trial);free(trans);return-1;}
            *out++=gs&0xFF;*out++=(gs>>8)&0xFF;*out++=(gs>>16)&0xFF;*out++=gs>>24;
            memcpy(out,gp,gs);out+=gs;
        }
    }

    h->magic=CPENC_MAGIC;h->version=CPENC_VERSION;h->stride=(uint8_t)stride;
    /* v2: store the stride remainder (tail bytes not covered by any group)
     * so decompress can restore them exactly */
    {
        size_t gend = gs * (size_t)stride;
        if (gend < sl) {
            size_t rem = sl - gend;
            if ((size_t)(edd - out) < rem) { free(trial); free(trans); return -1; }
            memcpy(out, trans + gend, rem);
            out += rem;
        }
    }
    h->original_size=(uint32_t)sl;h->compressed_size=(uint32_t)(out-(d+CPENC_HEADER_SIZE));
    *dl=(size_t)(out-d); free(trial); free(trans); return 0;
}

int cpenc_decompress(const void *sv, size_t sl, void *dv, size_t *dl) {
    const uint8_t *s=(const uint8_t*)sv; uint8_t *d=(uint8_t*)dv;
    const CPEncHeader *h=(const CPEncHeader*)s;
    if(sl<CPENC_HEADER_SIZE) return -1;
    if(h->magic!=CPENC_MAGIC||h->version!=CPENC_VERSION) return -1;
    if(h->stride==CPENC_STRIDE_RAW){
        if(h->original_size==0){*dl=0;return 0;}
        if(*dl<h->original_size) return -1;
        memcpy(d,s+CPENC_HEADER_SIZE,h->original_size);*dl=h->original_size;return 0;}

    int stride=(int)h->stride; if(stride<1||stride>8) return -1;
    size_t orig=h->original_size, gs=orig/(size_t)stride;
    if(*dl<orig) return -1;
    const uint8_t *p=s+CPENC_HEADER_SIZE, *edd=s+sl;
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
        }else{
            if(h_dec(&p,edd,gp,gs)<0){if(stride>1)free(trans);return-1;}
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
    *dl=orig; return 0;
}
