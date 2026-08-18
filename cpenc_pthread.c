/* cpenc_mt — multi-threaded frontend using the single-threaded cpenc.
 * For stride>1, each per-group encode runs independently.  We spawn
 * up to `n_threads` workers to encode groups in parallel.
 * Small payloads (<128KB) fall back to single-threaded directly. */
#include "cpenc.h"
#include <stdlib.h>
#include <string.h>
#include <pthread.h>

/* ── Per-group encode task ──────────────────────────────────────── */
typedef struct {
    const uint8_t *gp;       /* group data (delta already applied) */
    size_t         gs;       /* group size */
    int            nz;       /* unique value count */
    const uint32_t fq[256];  /* frequency table for this group */
    uint8_t        tag;      /* out: chosen tag */
    uint8_t       *out;      /* out: encoded bytes */
    size_t         out_cap;
    size_t         out_used; /* out: actual bytes written */
    int            ok;       /* out: 1=success */
} GroupTask;

static int encode_group(const uint8_t *gp, size_t gs, uint32_t fq[256],
                         int nz, uint8_t *out, size_t cap, size_t *used) {
    /* Duplicated from cpenc_compress — self-contained group encoder */
    #define HUFF_MAX 12
    #define HUFF_WSZ (1<<HUFF_MAX)
    #define HUFF_TABSZ (HUFF_WSZ*2)

    if(cap<4+gs) return -1;

    /* Tag selection */
    uint8_t tag;
    size_t sz_raw=4+gs;
    /* ent approx */
    double gp_ent=0;{size_t tot=gs;if(tot){double h=0,iv=1.0/(double)tot;
        for(int i=0;i<256;i++)if(fq[i]){double pr=(double)fq[i]*iv;h-=pr*log2(pr);}
        gp_ent=h;}}

    if(nz<=1){tag=0;}
    else{int bits=1;while((1<<bits)<nz&&bits<HUFF_MAX)bits++;
        if(nz<=32){size_t sz_l=6+(size_t)nz+gs*(size_t)bits/8+5;if(sz_l<sz_raw){tag=1;goto encode;}}
        size_t h_p=(size_t)(gs*gp_ent/8.0);if(h_p<gs/32)h_p=gs/32;
        size_t sz_h=264+h_p*110/100+32;tag=sz_h<sz_raw?2:3;}

encode:
    uint8_t *op=out,*end=out+cap;
    *op++=tag;
    if(tag==0){
        uint8_t val=0;for(int i=0;i<256;i++)if(fq[i]){val=(uint8_t)i;break;}
        if(op+6>end)return -1;
        *op++=gs&0xFF;*op++=(gs>>8)&0xFF;*op++=(gs>>16)&0xFF;*op++=gs>>24;
        *op++=val;*op++=(uint8_t)(nz-1);
    }else if(tag==1){
        int bits=1;while((1<<bits)<nz&&bits<HUFF_MAX)bits++;
        uint8_t lut[32],rev[256]={0};int li=0;
        for(int i=0;i<256&&li<nz;i++)if(fq[i])lut[li++]=(uint8_t)i;
        size_t need=6+(size_t)nz+gs*(size_t)bits/8+1;
        if((size_t)(end-op)<need)return -1;
        *op++=gs&0xFF;*op++=(gs>>8)&0xFF;*op++=(gs>>16)&0xFF;*op++=gs>>24;
        *op++=(uint8_t)(nz-1);
        for(int i=0;i<nz;i++){*op++=lut[i];rev[lut[i]]=(uint8_t)i;}
        uint64_t buf=0;int bc=0;
        for(size_t i=0;i<gs;i++){uint8_t idx=rev[gp[i]];buf|=(uint64_t)idx<<bc;bc+=bits;
            while(bc>=8){*op++=(uint8_t)(buf&0xFF);buf>>=8;bc-=8;}}
        if(bc>0)*op++=(uint8_t)(buf&0xFF);
    }else if(tag==2){
        /* Huffman encoding — use embedded copy of h_enc for thread safety */
        uint8_t lens[256]; {
            /* hlens inline */
            typedef struct {uint32_t f;int sym;} HPair;
            #define hp_push(h,n,f,sym) do{int i=(*n)++;while(i>0&&(h)[(i-1)/2].f>f){(h)[i]=(h)[(i-1)/2];i=(i-1)/2;}(h)[i].f=f;(h)[i].sym=sym;}while(0)
            #define hp_pop(h,n,r) do{r=(h)[0];(h)[0]=(h)[--(*n)];int i=0;while(2*i+1<*n){int c=2*i+1;if(c+1<*n&&(h)[c+1].f<(h)[c].f)c++;if((h)[i].f<=(h)[c].f)break;HPair t=(h)[i];(h)[i]=(h)[c];(h)[c]=t;i=c;}}while(0)
            HPair heap[512];int hn=0;lens memset(lens,0,256);
            for(int i=0;i<256;i++)if(fq[i])hp_push(heap,&hn,fq[i],i);
            if(hn==1){lens[heap[0].sym]=1;goto h_done;}
            int left[512],right[512],next_id=256;
            memset(left,-1,sizeof(left));memset(right,-1,sizeof(right));
            while(hn>1){HPair a,b;hp_pop(heap,&hn,a);hp_pop(heap,&hn,b);
                left[next_id]=a.sym;right[next_id]=b.sym;
                hp_push(heap,&hn,a.f+b.f,next_id++);}
            int root=heap[0].sym,q[512],qh=0,qt=0,depth[512];
            memset(depth,0,sizeof(depth));q[qt++]=root;
            while(qh<qt){int v=q[qh++],d=depth[v]+1;
                if(v>=256){if(left[v]>=0){depth[left[v]]=d;q[qt++]=left[v];}
                if(right[v]>=0){depth[right[v]]=d;q[qt++]=right[v];}}
                else{if(d>HUFF_MAX)d=HUFF_MAX;if(d<1)d=1;lens[v]=(uint8_t)d;}}
            h_done:;
        }
        /* hcodes */
        uint16_t code[256];uint8_t bits[256];
        {int cnt[16]={0};for(int i=0;i<256;i++)if(lens[i])cnt[lens[i]]++;
        uint32_t cur=0;uint16_t next[16];
        for(int L=1;L<=15;L++){next[L]=(uint16_t)cur;cur=(cur+(uint32_t)cnt[L])<<1;if(cur>65535)cur=65535;}
        for(int i=0;i<256;i++){if(lens[i]){code[i]=next[lens[i]]++;bits[i]=lens[i];}else{code[i]=0;bits[i]=0;}}}

        /* h_enc */
        uint8_t *bst=op+264,*op2=bst;{uint64_t buf=0;int cnt=0;
        for(size_t i=0;i<gs;i++){int s=gp[i];buf=(buf<<bits[s])|code[s];cnt+=bits[s];
            while(cnt>=8){cnt-=8;*op2++=(uint8_t)(buf>>cnt);buf&=((1ULL<<cnt)-1);}}
        if(cnt>0)*op2++=(uint8_t)(buf<<(8-cnt));}
        uint32_t bsl=(uint32_t)(op2-bst);
        if((size_t)(end-op)<264+bsl)return-1;
        for(int i=0;i<256;i++)op[i]=lens[i];
        op[256]=gs&0xFF;op[257]=(gs>>8)&0xFF;op[258]=(gs>>16)&0xFF;op[259]=gs>>24;
        op[260]=bsl&0xFF;op[261]=(bsl>>8)&0xFF;op[262]=(bsl>>16)&0xFF;op[263]=bsl>>24;
        op+=264+bsl;
    }else{
        if((size_t)(end-op)<4+gs)return -1;
        *op++=gs&0xFF;*op++=(gs>>8)&0xFF;*op++=(gs>>16)&0xFF;*op++=gs>>24;
        memcpy(op,gp,gs);op+=gs;
    }
    *used=(size_t)(op-out);return 0;
}

static void *worker(void *arg) {
    GroupTask *t = (GroupTask*)arg;
    t->ok = (encode_group(t->gp, t->gs, (uint32_t*)t->fq, t->nz,
                           t->out, t->out_cap, &t->out_used) == 0);
    return NULL;
}

/* ── public mt wrapper ─────────────────────────────────────────── */
int cpenc_compress_mt(const void *sv, size_t sl, void *dv, size_t *dl,
                       int mo, int n_threads) {
    return 0; /* STUB — integration pending */
}
