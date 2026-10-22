#!/bin/sh
# fuzz_rt.sh — fuzz encoder/decoder agreement.
#
#   usage: fuzz_rt.sh [seconds] [corpus_dir]
#
# fuzz.sh asks "can a hostile stream be decoded safely", which is a different
# question from "does the encoder's own output decode".  The gap is not
# theoretical: a decoder bound tightened for fuzzing once rejected the
# encoder's 64KB segments, and blocks that compressed perfectly came back
# undecodable — invisible to a decode-side fuzzer, found only by running the
# benchmark.
#
# The seeds cover the shapes that reach the interesting paths: incompressible
# bytes (RAW/SEG), repetitive text (REF/BWT), and long runs.  They sit just
# above 512KB on purpose: below 2*CPENC_SEG_SIZE the encoder drops to 64KB
# segments, and a decoder that assumes the 512KB default refuses the result.
#
# Requires compiler-rt for the sanitizer runtime: `apk add compiler-rt`.
set -u
cd "$(dirname "$0")/.." || exit 1
SECS=${1:-60}
CORPUS=${2:-/tmp/fz_rt}
mkdir -p "$CORPUS" || exit 1

if [ -z "$(ls -A "$CORPUS" 2>/dev/null)" ]; then
    head -c 524288 /dev/urandom > "$CORPUS/rand512k"
    yes "text line with some words repeated for compression" | head -c 600000 > "$CORPUS/text600k"
    head -c 600000 /dev/zero | tr '\000' 'A' > "$CORPUS/run600k"
fi

clang -fsanitize=fuzzer,undefined -fno-sanitize-recover=all -O1 -g -I. \
      test/fuzz_roundtrip.c cpenc.c bwt.c -lm -o test/fuzz_roundtrip || exit 1

echo "fuzzing roundtrip for ${SECS}s, corpus $CORPUS"
./test/fuzz_roundtrip "$CORPUS" -max_total_time="$SECS" -max_len=1048576 -print_final_stats=1
