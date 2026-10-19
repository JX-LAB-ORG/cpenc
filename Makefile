.PHONY: test clean all

# Target-specific ISA flags: NEON on aarch64, a conservative x86-64 baseline
# elsewhere (the NEON paths are compiled out there).  Override with ARCHFLAGS=.
ARCH := $(shell uname -m)
ifeq ($(ARCH),aarch64)
ARCHFLAGS ?= -march=armv8-a+simd
else
ARCHFLAGS ?= -march=x86-64-v2
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

test: cpenc.o cpool.o test/cpenc_test.c test/cpenc_roundtrip_test.c test/cpool_test.c test/cpool_disk_test.c
	cc $(CFLAGS) test/cpenc_test.c cpenc.o bwt.o -lm -o test/cpenc_test
	./test/cpenc_test
	cc $(CFLAGS) test/cpenc_roundtrip_test.c cpenc.o bwt.o -lm -o test/cpenc_roundtrip_test
	./test/cpenc_roundtrip_test
	cc $(CFLAGS) test/cpool_test.c cpool.o cpenc.o bwt.o -lm -lpthread -o test/cpool_test
	./test/cpool_test
	cc $(CFLAGS) test/cpool_disk_test.c cpool.o cpenc.o bwt.o -lm -lpthread -o test/cpool_disk_test
	./test/cpool_disk_test

clean:
	rm -f *.o *.so test/cpenc_test test/cpool_test
