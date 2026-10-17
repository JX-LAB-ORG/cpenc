#include "../cpenc.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int failed=0;
#define CHK(c,f,...) do{if(!(c)){fprintf(stderr,"  FAIL: "f"\n",##__VA_ARGS__);failed++;}}while(0)

static int rt(const char *nm, const uint8_t *d, size_t n, int mo) {
    size_t cap=cpenc_max_compressed_size(n);
    uint8_t *c=malloc(cap),*x=malloc(n+16);
    if(!c||!x){free(c);free(x);return-1;}
    size_t cl=cap,xl=n+16;
    int rc=cpenc_compress(d,n,c,&cl,mo);
    CHK(rc==0,"%s compress rc=%d",nm,rc);
    if(rc){free(c);free(x);return-1;}
    rc=cpenc_decompress(c,cl,x,&xl);
    CHK(rc==0,"%s decomp rc=%d",nm,rc);
    int ok=(rc==0&&xl==n&&memcmp(d,x,n)==0);
    CHK(ok,"%s roundtrip",nm);
    CPEncHeader*h=(CPEncHeader*)c;
    double ratio=n>0?(double)n/(double)cl:1;
    int raw=(h->stride==CPENC_STRIDE_RAW);
    fprintf(stderr,"  %s: %zu→%zuB (%.2fx) s=%u%s\n",nm,n,cl,ratio,(unsigned)h->stride,raw?" RAW":"");
    free(c);free(x);
    return ok?0:-1;
}

static void gl(uint8_t*b,size_t n){for(size_t i=0;i<n;i++)b[i]=(uint8_t)((i*1000003+17)&0xFF);}
static void gi(uint8_t*b,size_t n){for(size_t i=0;i<n;i++)b[i]=i&1?((uint8_t)((i/2)%16)<<4)|((uint8_t)((i/2+1)%16)):(uint8_t)((i/2)%16)|((uint8_t)((i/2+1)%16)<<4);}
static void gc(uint8_t*b,size_t n){memset(b,0x42,n);for(size_t i=0;i<n;i+=8)b[i]=0x42+(uint8_t)(i&3);}
static void gs(uint8_t*b,size_t n){for(size_t i=0;i<n;i++)b[i]=(i*7)&0xFF;}
static void gr(uint8_t*b,size_t n){for(size_t i=0;i<n;i++)b[i]=rand()&0xFF;}
static void gt(uint8_t*b,size_t n){const char*t="The quick brown fox jumps over the lazy dog. ";size_t tl=strlen(t);for(size_t i=0;i<n;i++)b[i]=(uint8_t)t[i%tl];}
/* Reference layer: exact repeats (identical lines, records with varying
 * fields, random data with embedded duplicate blocks) must roundtrip. */
static void gref1(uint8_t*b,size_t n){const char*l="2026-10-16 ERROR something failed with code 42\n";size_t L=strlen(l);for(size_t i=0;i<n;i++)b[i]=(uint8_t)l[i%L];}
static void gref2(uint8_t*b,size_t n){size_t i=0;unsigned v=0;while(i<n){char t[96];int L=snprintf(t,sizeof(t),"rec %08u id=%08x status=ok\n",v,v*2654435761u);v++;for(int k=0;k<L&&i<n;k++)b[i++]=(uint8_t)t[k];}}
static void gref3(uint8_t*b,size_t n){for(size_t i=0;i<n;i++)b[i]=(uint8_t)((i*2654435761u)>>24);
    for(size_t k=0;k<8;k++){size_t blk=200+k*40,src=(k*997)%(n-blk),dst=(k*1777)%(n-blk);memmove(b+dst,b+src,blk);}}
/* Approximate match ("base + exceptions"): long records differing in a few
 * scattered bytes.  These look like noise to delta/BWT (byte entropy ~8)
 * and used to fall to RAW at 1.00x. */
static void gref4(uint8_t*b,size_t n){size_t RL=1024; uint8_t*rec=malloc(RL); uint64_t s=12345;
    for(size_t i=0;i<RL;i++){s^=s<<13;s^=s>>7;s^=s<<17;rec[i]=(uint8_t)(s>>32);}
    for(size_t o=0;o+RL<=n;o+=RL){memcpy(b+o,rec,RL);
        for(int k=0;k<3;k++){s^=s<<13;s^=s>>7;s^=s<<17;size_t p=(s>>32)%RL;
            s^=s<<13;s^=s>>7;s^=s<<17;rec[p]=(uint8_t)(s>>32);}}
    free(rec);}
static void gref5(uint8_t*b,size_t n){size_t RL=4096; uint8_t*rec=malloc(RL); uint64_t s=999;
    for(size_t i=0;i<RL;i++){s^=s<<13;s^=s>>7;s^=s<<17;rec[i]=(uint8_t)(s>>32);}
    for(size_t o=0;o+RL<=n;o+=RL){memcpy(b+o,rec,RL);
        for(int k=0;k<8;k++){s^=s<<13;s^=s>>7;s^=s<<17;size_t p=(s>>32)%RL;
            s^=s<<13;s^=s>>7;s^=s<<17;rec[p]=(uint8_t)(s>>32);}}
    free(rec);}

int main(void){
    uint8_t*b=malloc(32768);if(!b)return 1;
    srand(42);

    fprintf(stderr,"=== CPEnc v0.1 ===\n");
    typedef void(*gf)(uint8_t*,size_t);
    struct{const char*nm;gf fn;int mo;}ts[]={
        {"INT4_4K",gi,2},{"LINEAR_4K",gl,2},{"CONST_4K",gc,1},
        {"SLOPE_4K",gs,2},{"TEXT_4K",gt,2},{"INT4_16K",gi,4},
        {"LINEAR_16K",gl,4},{"REFREP_16K",gref1,2},{"REFREC_16K",gref2,2},
        {"REFMIX_16K",gref3,2},{"REFREC1K_16K",gref4,2},
        {"REFREC4K_16K",gref5,2},{0,0,0}};
    for(int i=0;ts[i].nm;i++){
        size_t sz=strstr(ts[i].nm,"16K")?16384:4096;
        ts[i].fn(b,sz); rt(ts[i].nm,b,sz,ts[i].mo);
    }

    fprintf(stderr,"=== Size Sweep ===\n");
    size_t szs[]={64,256,1024,4096,16384,0};
    for(int i=0;szs[i];i++){gl(b,szs[i]);char nm[32];snprintf(nm,sizeof(nm),"L_%zu",szs[i]);rt(nm,b,szs[i],2);}

    fprintf(stderr,"=== RAW ===\n");
    gr(b,4096);
    {size_t cap=cpenc_max_compressed_size(4096);uint8_t*c=malloc(cap),*x=malloc(4112);size_t cl=cap,xl=4112;
     cpenc_compress(b,4096,c,&cl,2);CPEncHeader*h=(CPEncHeader*)c;
     fprintf(stderr,"  RND: s=%u c=%zu r=%.2fx\n",(unsigned)h->stride,cl,(double)4096/cl);
     cpenc_decompress(c,cl,x,&xl);CHK(xl==4096&&memcmp(b,x,4096)==0,"RAW");free(c);free(x);}

    fprintf(stderr,"=== Edge ===\n");
    {uint8_t e[256],x[256];size_t cl=sizeof(e),xl=sizeof(x);
     CHK(cpenc_compress(NULL,0,e,&cl,2)==0,"empty");
     CHK(cl>=CPENC_HEADER_SIZE,"empty hdr");
     CHK(cpenc_decompress(e,cl,x,&xl)==0,"empty dec");
     CHK(xl==0,"empty len");}
    {uint8_t e[256];gl(e,100);size_t cl=12;
     CHK(cpenc_compress(e,100,e,&cl,2)==-1,"toosmall");}
    {uint8_t e[64],x[64];memset(e,0xAA,64);size_t xl=64;
     CHK(cpenc_decompress(e,64,x,&xl)==-1,"badmagic");}

    fprintf(stderr,"\n=== %d failures ===\n",failed);
    free(b);
    return failed?1:0;
}
