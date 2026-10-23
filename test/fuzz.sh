#!/bin/sh
# fuzz.sh — fuzz the decoder with libFuzzer, under UBSan and guard pages.
#
#   usage: fuzz.sh [seconds] [corpus_dir]
#
# AddressSanitizer does not run in this environment (PRoot's address-space
# emulation rejects its shadow map), so the target places both buffers
# against PROT_NONE guard pages and UBSan covers undefined arithmetic.  See
# test/fuzz_cpenc.c.
#
# Requires compiler-rt for the sanitizer runtime: `apk add compiler-rt`.
set -u
cd "$(dirname "$0")/.." || exit 1
SECS=${1:-60}
CORPUS=${2:-/tmp/fz_cpenc}
mkdir -p "$CORPUS" || exit 1

cc -std=c11 -D_GNU_SOURCE -Wall -O2 -I. test/fuzz_seed.c cpenc.o bwt.o -lm -o test/fuzz_seed || exit 1
if [ -z "$(ls -A "$CORPUS" 2>/dev/null)" ]; then
    ./test/fuzz_seed "$CORPUS" || exit 1
fi

# Hand-built seeds that random mutation cannot construct (a self-consistent
# BWT block whose capacity prefix is deliberately too small) are kept in
# test/fuzz_seeds/ and copied into the corpus on every run, including runs
# against a corpus that is already populated.
if [ -d test/fuzz_seeds ]; then
    for f in test/fuzz_seeds/*.seed; do
        [ -e "$f" ] || continue
        cp "$f" "$CORPUS"/ || exit 1
    done
fi

clang -fsanitize=fuzzer,undefined -fno-sanitize-recover=all -O1 -g -I. \
      test/fuzz_cpenc.c cpenc.c bwt.c -lm -o test/fuzz_cpenc || exit 1

echo "fuzzing for ${SECS}s, corpus $CORPUS"
./test/fuzz_cpenc "$CORPUS" -max_total_time="$SECS" -max_len=131072 -print_final_stats=1
