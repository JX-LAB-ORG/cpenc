/* fuzz_cpenc.c — libFuzzer target for the block decoder.
 *
 * In the container, a block's bytes come off a disk, so every field in the
 * CPEnc header and every length inside the stream is untrusted.  The
 * contract cpenc_decompress() has to hold is narrow and testable: a
 * malformed block makes it return non-zero, it never writes past the *dl
 * capacity it was given, and it never reads past the sl bytes it was given.
 *
 * Why guard pages instead of ASan: AddressSanitizer does not run in this
 * environment (PRoot's address-space emulation rejects its shadow map), so
 * both buffers are placed flush against PROT_NONE pages.  An overrun traps
 * on the spot — the same signal ASan would produce, only at page
 * granularity.  UBSan is enabled alongside, so undefined arithmetic and
 * shifts are caught too.
 *
 * Input layout: [u32 capacity][block bytes].  The capacity is fuzzer
 * controlled on purpose — in the container it comes from the inode table,
 * which is just as untrusted as the payload, and the interesting bugs live
 * in the disagreement between the two.
 */
#include "../cpenc.h"
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <unistd.h>

#define PG 4096u

/* A buffer whose last byte sits immediately before an unmapped page, with
 * another unmapped page before it.  `len` is usable; the returned pointer is
 * where the caller's data starts. */
static uint8_t *guard_buf(size_t len, size_t *mapped_len, void **map_base) {
    size_t pages = (len + PG - 1) / PG;
    if (pages == 0) pages = 1;
    size_t total = (pages + 2) * PG;
    uint8_t *base = mmap(NULL, total, PROT_READ | PROT_WRITE,
                         MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    if (base == MAP_FAILED) return NULL;
    mprotect(base, PG, PROT_NONE);
    mprotect(base + (pages + 1) * PG, PG, PROT_NONE);
    *mapped_len = total;
    *map_base = base;
    return base + PG + pages * PG - len;   /* end flush against the guard */
}

int LLVMFuzzerTestOneInput(const uint8_t *d, size_t n) {
    if (n < 8 || n > (4u << 20)) return 0;

    uint32_t cap;
    memcpy(&cap, d, 4);
    cap &= 0x1FFFFFu;                     /* at most 2MB of output */
    const uint8_t *blk = d + 4;
    size_t blen = n - 4;

    /* input, guarded at both ends: a decoder that reads past blen traps */
    size_t inm, outm;
    void *inb, *outb;
    uint8_t *in = guard_buf(blen, &inm, &inb);
    uint8_t *out = guard_buf(cap ? cap : 1, &outm, &outb);
    if (!in || !out) { if (in) munmap(inb, inm); if (out) munmap(outb, outm); return 0; }
    memcpy(in, blk, blen);

    size_t dl = cap;
    cpenc_decompress(in, blen, out, &dl);

    /* A second pass with the capacity the header itself asks for: callers
     * that size their buffer from the same untrusted place have no slack,
     * and that is exactly the case where an unchecked write shows up. */
    if (blen >= CPENC_HEADER_SIZE) {
        CPEncHeader h;
        memcpy(&h, in, sizeof h);
        if (h.original_size && h.original_size <= (1u << 21)) {
            size_t dl2 = h.original_size;
            size_t m2; void *b2;
            uint8_t *o2 = guard_buf(dl2, &m2, &b2);
            if (o2) {
                cpenc_decompress(in, blen, o2, &dl2);
                munmap(b2, m2);
            }
        }
    }

    munmap(inb, inm);
    munmap(outb, outm);
    return 0;
}
