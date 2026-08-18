.PHONY: test clean all

CFLAGS = -std=c11 -Wall -O3 -march=armv8-a+simd -I.

all: cpenc.o cpool.o libcpool_hook.o

cpenc.o: cpenc.c cpenc.h
	cc $(CFLAGS) -c cpenc.c -o cpenc.o

cpool.o: cpool.c cpool.h cpenc.h
	cc $(CFLAGS) -c cpool.c -o cpool.o

libcpool_hook.o: libcpool_hook.c cpool.h
	cc -std=gnu11 $(CFLAGS) -fPIC -c libcpool_hook.c -o libcpool_hook.o

libcpool_hook.so: cpenc.o cpool.o libcpool_hook.o
	cc -std=gnu11 $(CFLAGS) -shared -fPIC -o $@ $^ -lpthread -ldl

test: cpenc.o cpool.o test/cpenc_test.c test/cpool_test.c test/cpool_disk_test.c
	cc $(CFLAGS) test/cpenc_test.c cpenc.o -lm -o test/cpenc_test
	./test/cpenc_test
	cc $(CFLAGS) test/cpool_test.c cpool.o cpenc.o -lm -lpthread -o test/cpool_test
	./test/cpool_test
	cc $(CFLAGS) test/cpool_disk_test.c cpool.o cpenc.o -lm -lpthread -o test/cpool_disk_test
	./test/cpool_disk_test

clean:
	rm -f *.o *.so test/cpenc_test test/cpool_test
