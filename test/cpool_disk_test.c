#include "../cpool.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

static int failed=0;
#define CHK(c,m) do{if(!(c)){failed++;fprintf(stderr,"  FAIL: %s\n",m);}else{fprintf(stderr,"  OK %s\n",m);}}while(0)

int main(void) {
    fprintf(stderr,"=== Disk Backing ===\n");
    /* Small warm RAM (64KB), large disk backing — forces spill */
    CPoolConfig cfg = {
        .warm_limit  = 65536,   /* 64KB warm */
        .file_path   = "/tmp/cpool_disk_test.dat",
        .file_limit  = 4ULL * 1024 * 1024,  /* 4MB disk */
    };
    unlink(cfg.file_path);
    CPool *p = cpool_create(&cfg);
    CHK(p!=NULL,"create with disk"); if(!p) return 1;

    /* Allocate blocks totalling > 64KB to force eviction to disk */
    CHandle handles[32];
    uint8_t patterns[32];
    for(int i=0;i<32;i++){
        handles[i]=cpool_alloc(p,4096);
        CHK(handles[i]!=0,"alloc 4K");
        patterns[i]=(uint8_t)((i*7+13)&0xFF);
        void *a=cpool_access(p,handles[i]);
        CHK(a!=NULL,"access");
        memset(a,patterns[i],4096);
        cpool_release(p,handles[i]);
    }

    /* Verify all 32 blocks survive roundtrip (some evicted to disk) */
    int ok=1;
    for(int i=0;i<32;i++){
        void *v=cpool_access(p,handles[i]);
        for(size_t j=0;j<4096;j++)
            if(((uint8_t*)v)[j]!=patterns[i]){ok=0;break;}
        cpool_release(p,handles[i]);
        if(!ok)break;
    }
    CHK(ok,"32-block disk roundtrip");

    /* Check stats */
    CPoolStats st;
    cpool_stats(p,&st);
    fprintf(stderr,"  warm=%zu disk=%zu file_bytes=%zu ratio=%.2f\n",
        st.warm_used, (size_t)st.file_write_bytes, st.file_write_bytes, st.compress_ratio);
    CHK(st.compress_count>=32,"compressed 32+ blocks");
    CHK(st.file_write_bytes>0||st.warm_used>0,"data persisted somewhere");

    /* Free all */
    for(int i=0;i<32;i++) cpool_free(p,handles[i]);

    /* Test: alloc large block that forces direct-to-disk */
    CHandle h=cpool_alloc(p,131072); /* 128KB */
    CHK(h!=0,"alloc 128K");
    void *b=cpool_access(p,h);
    CHK(b!=NULL,"access 128K");
    memset(b,0xA5,131072);
    cpool_release(p,h);
    void *c=cpool_access(p,h);
    int ok2=1;
    for(size_t j=0;j<131072;j++) if(((uint8_t*)c)[j]!=0xA5) ok2=0;
    CHK(ok2,"128K roundtrip (spill)");
    cpool_release(p,h);
    cpool_free(p,h);

    cpool_destroy(p);
    /* Verify disk file was cleaned up */
    CHK(access("/tmp/cpool_disk_test.dat",F_OK)!=0,"disk file cleaned");

    fprintf(stderr,"\n=== %d failures ===\n",failed);
    return failed?1:0;
}
