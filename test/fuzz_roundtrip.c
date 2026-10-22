/* fuzz_roundtrip.c — encoder/decoder agreement under the fuzzer.
 *
 * fuzz_cpenc checks that a hostile stream is decoded safely.  That is a
 * different question from whether the encoder's own output decodes, and the
 * gap is not theoretical: a decoder bound added while hardening the fuzzer
 * capped segmented blocks at original_size/512KB segments, the encoder emits
 * 64KB segments below 1MB, and blocks that compressed perfectly came back
 * undecodable.  No decode-side fuzzer can see that, and the round-trip tests
 * only covered the shapes someone thought to write down.
 *
 * So: take arbitrary bytes, compress them the way a caller would, and insist
 * that decompressing reproduces them.  A mismatch aborts, which the fuzzer
 * reports as a finding.
 */
#include "../cpenc.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

int LLVMFuzzerTestOneInput(const uint8_t *d, size_t n) {
    if (n < 2) return 0;
    int order = d[0] % 5;
    const uint8_t *src = d + 1;
    size_t sl = n - 1;

    size_t cap = cpenc_max_compressed_size(sl);
    uint8_t *cb = malloc(cap);
    uint8_t *db = malloc(sl);
    if (!cb || !db) { free(cb); free(db); return 0; }

    size_t cl = cap;
    if (cpenc_compress(src, sl, cb, &cl, order) == 0) {
        size_t dl = sl;
        int rc = cpenc_decompress(cb, cl, db, &dl);
        if (rc != 0) {
            fprintf(stderr, "ROUNDTRIP: decompress failed (order=%d n=%zu cl=%zu rc=%d)\n",
                    order, sl, cl, rc);
            abort();
        }
        if (dl != sl || memcmp(src, db, sl) != 0) {
            size_t k = 0;
            while (k < sl && k < dl && src[k] == db[k]) k++;
            fprintf(stderr, "ROUNDTRIP: mismatch (order=%d n=%zu dl=%zu first@%zu)\n",
                    order, sl, dl, k);
            abort();
        }
        /* The stream has to be deterministic too: the same input compressed
         * twice is the same bytes, or nothing can be hashed or compared. */
        uint8_t *cb2 = malloc(cap);
        if (cb2) {
            size_t cl2 = cap;
            if (cpenc_compress(src, sl, cb2, &cl2, order) == 0 &&
                (cl2 != cl || memcmp(cb, cb2, cl) != 0)) {
                fprintf(stderr, "ROUNDTRIP: encoder not deterministic (order=%d n=%zu)\n",
                        order, sl);
                free(cb2); free(cb); free(db);
                abort();
            }
            free(cb2);
        }
    }

    free(cb); free(db);
    return 0;
}
