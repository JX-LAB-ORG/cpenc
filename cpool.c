#include "cpool.h"
#include "cpenc.h"
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <fcntl.h>
#include <unistd.h>

#define MAX_BLOCKS    4096
#define HANDLE_SHIFT   12
#define HANDLE_MASK    ((1ULL << HANDLE_SHIFT) - 1)

typedef struct {
    uint64_t id;
    size_t   original_size;
    uint8_t  tier;         /* 1=warm(RAM), 2=disk */
    uint8_t  access_count;
    uint64_t last_access;
    void    *warm_ptr;
    size_t   warm_size;
    size_t   warm_alloc;
    void    *temp_buf;
    uint64_t file_offset;
    size_t   file_size;    /* compressed bytes on disk */
} BlockMeta;

typedef struct FlNode { struct FlNode *next; size_t size; } FlNode;

struct CPool {
    CPoolConfig     cfg;
    pthread_mutex_t lock;
    uint8_t        *warm_buf;
    size_t          warm_capacity;
    FlNode         *free_head;
    uint64_t        tick;
    BlockMeta       blocks[MAX_BLOCKS];
    uint32_t        n_blocks;
    CPoolStats      st;
    int             disk_fd;       /* disk backing file (-1 = backend mode) */
    void           *backend_ctx;   /* custom backend context */
    uint64_t        disk_end;      /* next write position */
};

static inline size_t align8(size_t x) { return (x+7)&~7ULL; }

static void fl_insert(struct CPool *p, void *ptr, size_t sz) {
    if (sz < sizeof(FlNode)) return;
    FlNode *n = (FlNode*)ptr; n->size = sz;
    FlNode **pp = &p->free_head;
    while (*pp && (uintptr_t)(*pp) < (uintptr_t)n) pp = &(*pp)->next;
    n->next = *pp; *pp = n;
    if (n->next && (uint8_t*)n + n->size == (uint8_t*)n->next)
        { n->size += n->next->size; n->next = n->next->next; }
    for (FlNode *cur = p->free_head; cur && cur != n; cur = cur->next)
        if ((uint8_t*)cur + cur->size == (uint8_t*)n)
            { cur->size += n->size; cur->next = n->next; break; }
}

static void *fl_alloc(struct CPool *p, size_t need, size_t *out_size) {
    need = align8(need);
    FlNode **pp = &p->free_head;
    while (*pp) {
        if ((*pp)->size >= need) {
            FlNode *n = *pp; *pp = n->next;
            size_t rem = n->size - need;
            if (rem >= sizeof(FlNode)+8) { fl_insert(p, (uint8_t*)n+need, rem); *out_size = need; }
            else *out_size = n->size;
            return n;
        }
        pp = &(*pp)->next;
    }
    return NULL;
}

/* ── Write compressed block to disk backing file ──────────────── */
static int disk_write(struct CPool *p, void *data, size_t sz, uint64_t *out_off) {
    if (p->cfg.file_limit > 0 && p->disk_end + sz > p->cfg.file_limit) return -1;
    *out_off = p->disk_end;
    if (p->cfg.backend_open) {
        if (p->cfg.backend_write(p->backend_ctx, data, sz, p->disk_end) != 0) return -1;
    } else {
        if (p->disk_fd < 0) return -1;
        if (pwrite(p->disk_fd, data, sz, (off_t)p->disk_end) != (ssize_t)sz) return -1;
    }
    p->disk_end += sz;
    p->st.file_write_bytes += sz;
    return 0;
}

/* ── Read compressed block from disk backing file ─────────────── */
static void *disk_read(struct CPool *p, uint64_t off, size_t sz) {
    void *buf = malloc(sz);
    if (!buf) return NULL;
    if (p->cfg.backend_open) {
        if (p->cfg.backend_read(p->backend_ctx, buf, sz, off) != 0) { free(buf); return NULL; }
    } else {
        if (p->disk_fd < 0) { free(buf); return NULL; }
        if (pread(p->disk_fd, buf, sz, (off_t)off) != (ssize_t)sz) { free(buf); return NULL; }
    }
    p->st.file_read_bytes += sz;
    return buf;
}

/* ── Evict one LRU block from warm RAM (tier=1).  If disk backing
     is configured, spill compressed data to disk first. ──────── */
static int evict_one(struct CPool *p) {
    uint64_t oldest = UINT64_MAX; int idx = -1;
    for (int i=0; i<MAX_BLOCKS; i++) {
        BlockMeta *b=&p->blocks[i];
        if (b->tier==1 && b->access_count==0 && b->warm_ptr && b->last_access<oldest)
            {oldest=b->last_access; idx=i;}
    }
    if (idx<0) return 0;
    BlockMeta *b=&p->blocks[idx];

    /* Spill to disk if backing configured (never in RAM-only mode) */
    if (!(p->cfg.flags & CPOOL_NO_EVICT) && (p->cfg.backend_open || p->disk_fd >= 0)) {
        uint64_t off;
        if (disk_write(p, b->warm_ptr, b->warm_size, &off) == 0) {
            b->file_offset = off;
            b->file_size   = b->warm_size;
            b->tier        = 2; /* now disk-resident */
            p->st.file_uncompressed += b->original_size;
        }
    }

    p->st.warm_used          -= b->warm_size;
    p->st.warm_uncompressed  -= b->original_size;
    if (b->warm_ptr) fl_insert(p, b->warm_ptr, b->warm_alloc);
    b->warm_ptr  = NULL;
    b->warm_size = 0;
    return 1;
}

static int make_space(struct CPool *p, size_t need) {
    /* Try free list first */
    for (FlNode *n=p->free_head;n;n=n->next) if(n->size>=need) return 1;
    /* RAM-only mode: never evict */
    if (p->cfg.flags & CPOOL_NO_EVICT) return 0;
    /* Evict up to 8 blocks */
    for (int i=0;i<8;i++) {
        if(!evict_one(p)) break;
        for (FlNode *n=p->free_head;n;n=n->next) if(n->size>=need) return 1;
    }
    return 0;
}

/* ── Lifecycle ────────────────────────────────────────────── */
CPool *cpool_create(const CPoolConfig *cfg) {
    CPool *p=calloc(1,sizeof(*p)); if(!p) return NULL;
    p->cfg=*cfg;
    if(p->cfg.warm_limit==0) p->cfg.warm_limit=256ULL*1024*1024;
    p->warm_capacity=p->cfg.warm_limit;
    p->warm_buf=malloc(p->warm_capacity); if(!p->warm_buf){free(p);return NULL;}
    fl_insert(p,p->warm_buf,p->warm_capacity);

    /* Open disk backing (custom backend takes precedence over file;
     * RAM-only mode never touches disk) */
    p->disk_fd = -1;
    if (p->cfg.flags & CPOOL_NO_EVICT) {
        /* no backing at all */
    } else if (p->cfg.backend_open) {
        if (p->cfg.backend_open(&p->backend_ctx, p->cfg.file_limit) != 0) {
            free(p->warm_buf); free(p); return NULL;
        }
        p->disk_end = 0;
    } else if (p->cfg.file_path && p->cfg.file_path[0]) {
        p->disk_fd = open(p->cfg.file_path, O_CREAT|O_RDWR|O_TRUNC, 0600);
        if (p->disk_fd < 0) { free(p->warm_buf); free(p); return NULL; }
        /* Pre-allocate sparse backing to prevent ENOSPC surprises */
        if (p->cfg.file_limit > 0)
            ftruncate(p->disk_fd, (off_t)p->cfg.file_limit);
        p->disk_end = 0;
    }

    pthread_mutex_init(&p->lock,NULL);
    return p;
}

void cpool_destroy(CPool *p) {
    if(!p) return;
    pthread_mutex_destroy(&p->lock);
    free(p->warm_buf);
    if (p->cfg.backend_open) {
        if (p->cfg.backend_close) p->cfg.backend_close(p->backend_ctx);
    } else if (p->disk_fd >= 0) {
        close(p->disk_fd);
        if (p->cfg.file_path && p->cfg.file_path[0])
            unlink(p->cfg.file_path);
    }
    free(p);
}

/* ── Alloc ────────────────────────────────────────────────── */
CHandle cpool_alloc(CPool *p, size_t size) {
    if(!p||size==0) return 0;
    pthread_mutex_lock(&p->lock);
    int slot=-1;
    for(int i=0;i<MAX_BLOCKS;i++) if(p->blocks[i].id==0){slot=i;break;}
    if(slot<0){pthread_mutex_unlock(&p->lock);return 0;}
    BlockMeta *b=&p->blocks[slot];
    b->id=((uint64_t)(++p->tick)<<HANDLE_SHIFT)|(uint64_t)(uint32_t)slot;
    b->original_size=size; b->tier=1; b->access_count=0;
    b->last_access=0; b->warm_ptr=NULL; b->warm_size=0; b->temp_buf=NULL;
    b->file_offset=0; b->file_size=0;
    p->n_blocks++;
    CHandle h=b->id;
    pthread_mutex_unlock(&p->lock);
    return h;
}

/* ── Free ─────────────────────────────────────────────────── */
void cpool_free(CPool *p, CHandle h) {
    if(!p||h==0) return;
    pthread_mutex_lock(&p->lock);
    uint32_t idx=(uint32_t)(h&HANDLE_MASK);
    if(idx>=MAX_BLOCKS){pthread_mutex_unlock(&p->lock);return;}
    BlockMeta *b=&p->blocks[idx];
    if(b->id!=h){pthread_mutex_unlock(&p->lock);return;}
    if(b->tier==1&&b->warm_ptr){
        p->st.warm_used-=b->warm_size; p->st.warm_uncompressed-=b->original_size;
        fl_insert(p,b->warm_ptr,b->warm_alloc);
    }
    if(b->tier==2)
        p->st.file_uncompressed -= b->original_size;
    free(b->temp_buf);
    memset(b,0,sizeof(*b)); p->n_blocks--;
    pthread_mutex_unlock(&p->lock);
}

/* ── Access ───────────────────────────────────────────────── */
void *cpool_access(CPool *p, CHandle h) {
    if(!p||h==0) return NULL;
    pthread_mutex_lock(&p->lock);
    uint32_t idx=(uint32_t)(h&HANDLE_MASK);
    if(idx>=MAX_BLOCKS){pthread_mutex_unlock(&p->lock);return NULL;}
    BlockMeta *b=&p->blocks[idx];
    if(b->id!=h){pthread_mutex_unlock(&p->lock);return NULL;}
    b->access_count++; b->last_access=p->tick++;
    if(b->temp_buf){pthread_mutex_unlock(&p->lock);return b->temp_buf;}
    b->temp_buf=malloc(b->original_size);
    if(!b->temp_buf){b->access_count--;pthread_mutex_unlock(&p->lock);return NULL;}

    if(b->tier==1 && b->warm_ptr){
        /* Decompress (or raw-copy) from warm RAM */
        size_t ol=b->original_size;
        if (p->cfg.flags & CPOOL_NO_COMPRESS) {
            memcpy(b->temp_buf, b->warm_ptr, b->original_size);
        } else if(cpenc_decompress(b->warm_ptr,b->warm_size,b->temp_buf,&ol)||ol!=b->original_size)
            {free(b->temp_buf);b->temp_buf=NULL;b->access_count--;pthread_mutex_unlock(&p->lock);return NULL;}
        p->st.decompress_count++;
    } else if(b->tier==2 && b->file_size>0){
        /* Read from disk, decompress (or raw-copy) */
        void *disk_buf = disk_read(p, b->file_offset, b->file_size);
        if(!disk_buf){free(b->temp_buf);b->temp_buf=NULL;b->access_count--;pthread_mutex_unlock(&p->lock);return NULL;}
        size_t ol=b->original_size;
        if (p->cfg.flags & CPOOL_NO_COMPRESS) {
            memcpy(b->temp_buf, disk_buf, b->original_size);
        } else if(cpenc_decompress(disk_buf,b->file_size,b->temp_buf,&ol)||ol!=b->original_size)
            {free(disk_buf);free(b->temp_buf);b->temp_buf=NULL;b->access_count--;pthread_mutex_unlock(&p->lock);return NULL;}
        free(disk_buf);
        p->st.decompress_count++;
    } else {
        /* Fresh allocation — zero-fill */
        memset(b->temp_buf,0,b->original_size);
    }
    pthread_mutex_unlock(&p->lock);
    return b->temp_buf;
}

/* ── Release ──────────────────────────────────────────────── */
void cpool_release(CPool *p, CHandle h) {
    if(!p||h==0) return;
    pthread_mutex_lock(&p->lock);
    uint32_t idx=(uint32_t)(h&HANDLE_MASK);
    if(idx>=MAX_BLOCKS){pthread_mutex_unlock(&p->lock);return;}
    BlockMeta *b=&p->blocks[idx];
    if(b->id!=h||b->access_count==0){pthread_mutex_unlock(&p->lock);return;}
    b->access_count--;
    if(b->access_count>0||!b->temp_buf){pthread_mutex_unlock(&p->lock);return;}

    /* Compress (or raw-copy in NO_COMPRESS mode) */
    size_t cl;
    uint8_t *comp;
    if (p->cfg.flags & CPOOL_NO_COMPRESS) {
        cl = b->original_size;
        comp = malloc(cl);
        if(!comp){pthread_mutex_unlock(&p->lock);return;}
        memcpy(comp, b->temp_buf, cl);
    } else {
        size_t mc=cpenc_max_compressed_size(b->original_size);
        comp=malloc(mc); if(!comp){pthread_mutex_unlock(&p->lock);return;}
        cl=mc;
        if(cpenc_compress(b->temp_buf,b->original_size,comp,&cl,2))
            {free(comp);pthread_mutex_unlock(&p->lock);return;}
    }

    /* Free old warm storage */
    if(b->warm_ptr){p->st.warm_used-=b->warm_size;p->st.warm_uncompressed-=b->original_size;
        fl_insert(p,b->warm_ptr,b->warm_alloc);}

    /* Try to allocate warm space; if not, evict */
    if(!make_space(p,cl)){
        for(int i=0;i<32;i++){if (p->cfg.flags & CPOOL_NO_EVICT) break; evict_one(p);if(make_space(p,cl))break;}
        if(!make_space(p,cl)){
            /* Cannot fit in warm RAM — store to disk directly (never in
             * RAM-only mode: the block stays resident in temp_buf) */
            if(!(p->cfg.flags & CPOOL_NO_EVICT) && (p->cfg.backend_open || p->disk_fd>=0)){
                uint64_t off;
                if(disk_write(p,comp,cl,&off)==0){
                    b->file_offset=off;b->file_size=cl;b->tier=2;
                    p->st.file_uncompressed+=b->original_size;
                    p->st.compress_count++;
                    free(b->temp_buf);b->temp_buf=NULL;free(comp);
                    pthread_mutex_unlock(&p->lock);return;
                }
            }
            free(comp);pthread_mutex_unlock(&p->lock);return;
        }
    }
    size_t alloc_sz; void *wp=fl_alloc(p,cl,&alloc_sz);
    if(!wp){free(comp);pthread_mutex_unlock(&p->lock);return;}
    memcpy(wp,comp,cl); free(comp);
    b->warm_ptr=wp;b->warm_size=cl;b->warm_alloc=alloc_sz;b->tier=1;
    p->st.warm_used+=cl;p->st.warm_uncompressed+=b->original_size;p->st.compress_count++;
    free(b->temp_buf);b->temp_buf=NULL;
    pthread_mutex_unlock(&p->lock);
}

void cpool_hint_cold(CPool *p, CHandle h) {
    if(!p||h==0) return;
    pthread_mutex_lock(&p->lock);
    uint32_t idx=(uint32_t)(h&HANDLE_MASK);
    if(idx<MAX_BLOCKS&&p->blocks[idx].id==h) p->blocks[idx].last_access=0;
    pthread_mutex_unlock(&p->lock);
}

void cpool_stats(CPool *p, CPoolStats *out) {
    if(!p||!out) return;
    pthread_mutex_lock(&p->lock);
    *out=p->st;
    out->compress_ratio= (p->st.warm_used+p->st.file_write_bytes) > 0
        ?(float)(p->st.warm_uncompressed+p->st.file_uncompressed)
         /(float)(p->st.warm_used+p->st.file_write_bytes)
        :1.0f;
    pthread_mutex_unlock(&p->lock);
}

size_t cpool_warm_free(const CPool *p) {
    if(!p) return 0;
    size_t total=0;
    for(FlNode *n=p->free_head;n;n=n->next) total+=n->size;
    return total;
}
