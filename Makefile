.PHONY: test clean all fuzz

# Target-specific ISA flags: NEON on aarch64, a conservative x86-64 baseline
# elsewhere (the NEON paths are compiled out there).  Override with ARCHFLAGS=.
ARCH := $(shell uname -m)
ifeq ($(ARCH),aarch64)
ARCHFLAGS ?= -march=armv8-a+simd
AVX2FLAGS ?=
else
ARCHFLAGS ?= -march=x86-64-v2
# x86-64-v2 has no AVX2, so the AVX2 kernel is off by default; the check
# target turns it on explicitly.  Appended after $(CFLAGS) so it wins.
AVX2FLAGS ?= -mavx2
endif
CFLAGS = -std=c11 -D_GNU_SOURCE -fPIC -Wall -O3 $(ARCHFLAGS) -I.

all: cpenc.o bwt.o cpool.o libcpool_hook.o

cpenc.o: cpenc.c cpenc.h bwt.h
	cc $(CFLAGS) -c cpenc.c -o cpenc.o

bwt.o: bwt.c bwt.h
	cc $(CFLAGS) -c bwt.c -o bwt.o

cpool.o: cpool.c cpool.h cpenc.h bwt.h
	cc $(CFLAGS) -c cpool.c -o cpool.o

libcpool_hook.o: libcpool_hook.c cpool.h
	cc $(CFLAGS) -std=gnu11 -fPIC -c libcpool_hook.c -o libcpool_hook.o

libcpool_hook.so: cpenc.o cpool.o libcpool_hook.o
	cc $(CFLAGS) -std=gnu11 -shared -fPIC -o $@ $^ -lpthread -ldl

test: cpenc.o cpool.o bwt.o test/cpenc_test.c test/cpenc_roundtrip_test.c test/cpool_test.c test/cpool_disk_test.c
	cc $(CFLAGS) test/cpenc_test.c cpenc.o bwt.o -lm -o test/cpenc_test
	./test/cpenc_test
	cc $(CFLAGS) test/cpenc_roundtrip_test.c cpenc.o bwt.o -lm -o test/cpenc_roundtrip_test
	./test/cpenc_roundtrip_test
	cc $(CFLAGS) test/cpool_test.c cpool.o cpenc.o bwt.o -lm -lpthread -o test/cpool_test
	./test/cpool_test
	cc $(CFLAGS) test/cpool_disk_test.c cpool.o cpenc.o bwt.o -lm -lpthread -o test/cpool_disk_test
	./test/cpool_disk_test
	cc $(CFLAGS) test/cpenc_rle_test.c cpenc.c bwt.o -lm -o test/cpenc_rle_test
	./test/cpenc_rle_test

# Codec benchmark: ratio + throughput over a real file, per stride.
bench: test/bench_cpenc test/cpenc_verify test/cpenc_bisect
test/bench_cpenc: test/bench_cpenc.c cpenc.o bwt.o
	cc $(CFLAGS) test/bench_cpenc.c cpenc.o bwt.o -lm -o test/bench_cpenc

# Round-trip sweep by block size; prints the failing offset and stride.
test/cpenc_verify: test/cpenc_verify.c cpenc.o bwt.o
	cc $(CFLAGS) test/cpenc_verify.c cpenc.o bwt.o -lm -o test/cpenc_verify

# Shrink a failing block to the smallest length that still fails.
test/cpenc_bisect: test/cpenc_bisect.c cpenc.o bwt.o
	cc $(CFLAGS) test/cpenc_bisect.c cpenc.o bwt.o -lm -o test/cpenc_bisect

# AVX2 LUT decode kernel vs the scalar reference.  Needs -mavx2 on x86-64;
# on aarch64 it builds the stub main (nothing to check there).
test/avx2_lut_check: test/avx2_lut_check.c cpenc_avx2.h
	cc $(CFLAGS) $(AVX2FLAGS) test/avx2_lut_check.c -o test/avx2_lut_check
	./test/avx2_lut_check

# Fuzz the decoder: libFuzzer + UBSan + guard pages (needs compiler-rt).
fuzz:
	sh test/fuzz.sh $(SECS)

# Fuzz encoder/decoder agreement: compress, decompress, insist on equality.
fuzz-rt:
	sh test/fuzz_rt.sh $(SECS)

clean:
	rm -f *.o *.so test/cpenc_test test/cpool_test test/fuzz_seed test/fuzz_cpenc test/fuzz_roundtrip test/bench_cpenc test/cpenc_verify test/cpenc_bisect test/cpenc_rle_test test/avx2_lut_check
