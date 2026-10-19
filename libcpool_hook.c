#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include "cpool.h"
#include <pthread.h>
#include <dlfcn.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>
#include <sys/mman.h>

/* ─── Threshold: only compress allocations >= 4KB to avoid overhead ── */
#define HOOK_THRESHOLD 4096

/* ─── Reentrancy guard: prevent infinite recursion when cpool ────────
       internally calls malloc/free (e.g. during cpool_create/access) */
static __thread int hook_busy;

/* ─── Global pool instance ──────────────────────────────────────────── */
static CPool *glob_pool;
static pthread_mutex_t glob_lock = PTHREAD_MUTEX_INITIALIZER;

/* ─── Pointer → CHandle mapping ──────────────────────────────────────── */
#define MAX_MAP 16384
typedef struct { void *ptr; CHandle h; size_t size; } PtrMap;
static PtrMap ptr_map[MAX_MAP];
static int ptr_map_n;

/* ─── Real function pointers ─────────────────────────────────────────── */
static void *(*real_malloc)(size_t)           = NULL;
static void  (*real_free)(void*)              = NULL;
static void *(*real_realloc)(void*, size_t)   = NULL;
static void *(*real_mmap)(void*, size_t, int, int, int, off_t) = NULL;
static int   (*real_munmap)(void*, size_t)    = NULL;

/* ─── Map helpers ────────────────────────────────────────────────────── */
static void map_add(void *ptr, CHandle h, size_t sz) {
    if (ptr_map_n >= MAX_MAP) return;
    for (int i = 0; i < MAX_MAP; i++) {
        if (ptr_map[i].ptr == NULL) {
            ptr_map[i].ptr  = ptr;
            ptr_map[i].h    = h;
            ptr_map[i].size = sz;
            ptr_map_n++;
            return;
        }
    }
}
static void map_del(void *ptr) {
    for (int i = 0; i < MAX_MAP; i++) {
        if (ptr_map[i].ptr == ptr) {
            ptr_map[i].ptr = NULL;
            ptr_map[i].h   = 0;
            ptr_map[i].size = 0;
            ptr_map_n--;
            return;
        }
    }
}
static CHandle map_find(void *ptr, size_t *sz) {
    for (int i = 0; i < MAX_MAP; i++) {
        if (ptr_map[i].ptr == ptr) {
            if (sz) *sz = ptr_map[i].size;
            return ptr_map[i].h;
        }
    }
    return 0;
}

/* ─── Lazy pool init ─────────────────────────────────────────────────── */
static CPool *get_pool(void) {
    if (!glob_pool) {
        CPoolConfig cfg = {
            .warm_limit  = 256ULL * 1024 * 1024, /* 256MB */
            .file_path   = NULL,
            .file_limit  = 0,
        };
        glob_pool = cpool_create(&cfg);
    }
    return glob_pool;
}

/* ═══════════════════════════════════════════════════════════════════════
   Hooked functions
   ═══════════════════════════════════════════════════════════════════════ */

void *malloc(size_t size) {
    if (!real_malloc)
        real_malloc = (typeof(real_malloc))dlsym(RTLD_NEXT, "malloc");
    if (hook_busy || size < HOOK_THRESHOLD)
        return real_malloc(size);

    hook_busy = 1;
    CPool *p = get_pool();
    if (!p) { hook_busy = 0; return real_malloc(size); }

    pthread_mutex_lock(&glob_lock);
    CHandle h = cpool_alloc(p, size);
    if (!h) { pthread_mutex_unlock(&glob_lock); hook_busy = 0; return real_malloc(size); }

    void *ptr = cpool_access(p, h);
    if (!ptr) { cpool_free(p, h); pthread_mutex_unlock(&glob_lock); hook_busy = 0; return real_malloc(size); }

    map_add(ptr, h, size);
    pthread_mutex_unlock(&glob_lock);
    hook_busy = 0;
    return ptr;
}

void free(void *ptr) {
    if (!ptr) return;
    if (!real_free)
        real_free = (typeof(real_free))dlsym(RTLD_NEXT, "free");
    if (hook_busy) { real_free(ptr); return; }

    hook_busy = 1;
    pthread_mutex_lock(&glob_lock);
    CHandle h = map_find(ptr, NULL);
    if (h && glob_pool) {
        map_del(ptr);
        cpool_release(glob_pool, h);
        cpool_free(glob_pool, h);
        pthread_mutex_unlock(&glob_lock);
        hook_busy = 0;
        return;
    }
    pthread_mutex_unlock(&glob_lock);
    hook_busy = 0;
    real_free(ptr);
}

void *realloc(void *ptr, size_t size) {
    if (!ptr) return malloc(size);
    if (size == 0) { free(ptr); return NULL; }
    if (!real_realloc)
        real_realloc = (typeof(real_realloc))dlsym(RTLD_NEXT, "realloc");
    if (hook_busy) return real_realloc(ptr, size);
    hook_busy = 1;

    pthread_mutex_lock(&glob_lock);
    size_t osz;
    CHandle h = map_find(ptr, &osz);
    if (h && glob_pool) {
        map_del(ptr);
        cpool_release(glob_pool, h);
        cpool_free(glob_pool, h);

        if (size < HOOK_THRESHOLD) {
            void *np = real_malloc(size);
            if (np && osz > 0) {
                size_t cp = size < osz ? size : osz;
                memcpy(np, ptr, cp);
            }
            pthread_mutex_unlock(&glob_lock);
            hook_busy = 0;
            return np;
        }

        CPool *p = get_pool();
        CHandle nh = cpool_alloc(p, size);
        if (!nh) { pthread_mutex_unlock(&glob_lock); hook_busy = 0; return real_realloc(ptr, size); }

        void *np = cpool_access(p, nh);
        if (!np) { cpool_free(p, nh); pthread_mutex_unlock(&glob_lock); hook_busy = 0; return real_realloc(ptr, size); }

        size_t cp = size < osz ? size : osz;
        if (cp > 0) memcpy(np, ptr, cp);

        map_add(np, nh, size);
        pthread_mutex_unlock(&glob_lock);
        hook_busy = 0;
        return np;
    }
    pthread_mutex_unlock(&glob_lock);
    hook_busy = 0;
    return real_realloc(ptr, size);
}

void *calloc(size_t nmemb, size_t size) {
    size_t total = nmemb * size;
    void *ptr = malloc(total);
    if (ptr) memset(ptr, 0, total);
    return ptr;
}

void *mmap(void *addr, size_t length, int prot, int flags,
           int fd, off_t offset) {
    if (!real_mmap)
        real_mmap = (typeof(real_mmap))dlsym(RTLD_NEXT, "mmap");
    if (hook_busy) return real_mmap(addr, length, prot, flags, fd, offset);
    hook_busy = 1;

    /* Anonymous mappings: compress via cpool */
    if (fd == -1 && (flags & MAP_ANONYMOUS) && length >= HOOK_THRESHOLD && glob_pool) {
        CPool *p = get_pool();
        if (!p) { hook_busy = 0; return (void*)-1 /* MAP_FAILED */; }

        pthread_mutex_lock(&glob_lock);
        CHandle h = cpool_alloc(p, length);
        if (!h) { pthread_mutex_unlock(&glob_lock); hook_busy = 0; return (void*)-1; }

        void *ptr = cpool_access(p, h);
        if (!ptr) {
            cpool_free(p, h);
            pthread_mutex_unlock(&glob_lock);
            hook_busy = 0;
            return (void*)-1;
        }
        map_add(ptr, h, length);
        pthread_mutex_unlock(&glob_lock);
        hook_busy = 0;
        return ptr;
    }
    hook_busy = 0;
    return real_mmap(addr, length, prot, flags, fd, offset);
}

int munmap(void *addr, size_t length) {
    if (!addr) return 0;
    if (!real_munmap)
        real_munmap = (typeof(real_munmap))dlsym(RTLD_NEXT, "munmap");
    if (hook_busy) return real_munmap(addr, length);
    hook_busy = 1;

    pthread_mutex_lock(&glob_lock);
    CHandle h = map_find(addr, NULL);
    if (h && glob_pool) {
        map_del(addr);
        cpool_release(glob_pool, h);
        cpool_free(glob_pool, h);
        pthread_mutex_unlock(&glob_lock);
        hook_busy = 0;
        return 0;
    }
    pthread_mutex_unlock(&glob_lock);
    hook_busy = 0;
    return real_munmap(addr, length);
}

/* ─── Environment-controlled pool size ────────────────────────────────── */
__attribute__((constructor)) static void _hook_init(void) {
    const char *e = getenv("CPOOL_LIMIT_MB");
    if (!e) e = getenv("LIBCPOOL_MB");
    if (e) {
        size_t mb = (size_t)atol(e);
        if (mb < 16) mb = 16;
        CPoolConfig cfg = { .warm_limit = mb * 1024 * 1024 };
        glob_pool = cpool_create(&cfg);
    }
    /* Otherwise lazy-init on first alloc */
}

__attribute__((destructor)) static void _hook_fini(void) {
    if (glob_pool) {
        cpool_destroy(glob_pool);
        glob_pool = NULL;
    }
}
