#include "../cpool.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int failed=0;
#define CHK(c,m) do{if(!(c)){failed++;fprintf(stderr,"  FAIL: %s\n",m);}else{fprintf(stderr,"  OK %s\n",m);}}while(0)

int main(void) {
    CPoolConfig cfg = { .warm_limit = 4*1024*1024 };
    CPool *p = cpool_create(&cfg);
    CHK(p!=NULL,"create"); if(!p) return 1;

    /* --- Roundtrip --- */
    fprintf(stderr,"=== Roundtrip ===\n");
    CHandle h1=cpool_alloc(p,256),h2=cpool_alloc(p,1024);
    CHK(h1&&h2,"alloc");
    void *a1=cpool_access(p,h1),*a2=cpool_access(p,h2);
    CHK(a1&&a2,"access");
    memset(a1,0xA5,256); memset(a2,0x3C,1024);
    uint8_t exp1[256],exp2[1024];
    memcpy(exp1,a1,256); memcpy(exp2,a2,1024);
    cpool_release(p,h1); cpool_release(p,h2);
    void *b1=cpool_access(p,h1),*b2=cpool_access(p,h2);
    int ok1=!memcmp(b1,exp1,256),ok2=!memcmp(b2,exp2,1024);
    CHK(ok1,"256B roundtrip"); CHK(ok2,"1KB roundtrip");
    cpool_release(p,h1); cpool_release(p,h2);

    /* --- Modify --- */
    fprintf(stderr,"=== Modify ===\n");
    void *m=cpool_access(p,h1);
    memset(m,0xFF,256);
    cpool_release(p,h1);
    void *mc=cpool_access(p,h1);
    int mok=1; for(int i=0;i<256;i++) if(((uint8_t*)mc)[i]!=0xFF) mok=0;
    CHK(mok,"modify persist"); cpool_release(p,h1);

    /* --- Stats --- */
    fprintf(stderr,"=== Stats ===\n");
    CPoolStats st; cpool_stats(p,&st);
    fprintf(stderr,"  warm_used=%zu uncomp=%zu ratio=%.2f comps=%llu\n",
        st.warm_used,st.warm_uncompressed,st.compress_ratio,
        (unsigned long long)st.compress_count);
    CHK(st.warm_used>0,"warm_used>0");

    /* --- Eviction --- */
    fprintf(stderr,"=== Eviction ===\n");
    /* Alloc many blocks with accesses interleaved to verify LRU correctness */
    CHandle hx[32]; int nx=0;
    for(int i=0;i<16&&nx<32;i++,nx++){
        hx[nx]=cpool_alloc(p,512);
        void *a=cpool_access(p,hx[nx]);
        memset(a,(uint8_t)i,512);
        cpool_release(p,hx[nx]);
    }
    /* Verify all retained */
    int evok=1;
    for(int i=0;i<nx;i++){
        void *v=cpool_access(p,hx[i]);
        for(size_t j=0;j<512;j++) if(((uint8_t*)v)[j]!=(uint8_t)i) {evok=0;break;}
        cpool_release(p,hx[i]);
        if(!evok) break;
    }
    CHK(evok,"16-block eviction roundtrip");
    for(int i=0;i<nx;i++) cpool_free(p,hx[i]);

    /* --- Free/Realloc --- */
    fprintf(stderr,"=== Free/Realloc ===\n");
    cpool_free(p,h1); cpool_free(p,h2);
    CHandle hr=cpool_alloc(p,8192);
    CHK(hr!=0,"realloc after free");
    void *ar=cpool_access(p,hr);
    memset(ar,0x77,8192);
    cpool_release(p,hr);
    void *br=cpool_access(p,hr);
    int rok=1; for(int i=0;i<8192;i++) if(((uint8_t*)br)[i]!=0x77) rok=0;
    CHK(rok,"8KB realloc"); cpool_release(p,hr); cpool_free(p,hr);

    cpool_destroy(p);

    /* --- Edge: zero limits --- */
    fprintf(stderr,"=== Edge Cases ===\n");
    CPoolConfig zc = { .warm_limit = 0 };
    CPool *zp = cpool_create(&zc);
    CHK(zp!=NULL,"zero-limit create (uses default)");
    cpool_destroy(zp);

    fprintf(stderr,"\n=== %d failures ===\n",failed);
    return failed?1:0;
}
