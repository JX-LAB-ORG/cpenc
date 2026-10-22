# CPEnc — the block codec behind the Container File System

Compress one block of data into a self-contained CPEnc stream, and get it back
byte for byte. No container, no filesystem, no daemon: memory in, memory out.

[libcfs](../libcfs), the Container File System (CFS), uses it to give 8MB
blocks random-access compression; it also works as a standalone library.

## What it is / what it is not

**It is** a compressor built for data that is structured but does not look like
text — numeric tables, binary records, log dumps, weight tensors. Per block it
picks one of several models (transpose, delta, BWT, references) and falls back
to RAW when the pick does not pay, so the worst case is **no compression**,
never **corruption**.

**It is not** a general-purpose compressor. On the same 673MB ext4 image:

| Tool | Ratio | Compress | Decompress |
|---|---:|---:|---:|
| **cpenc** | **3.68x** | **4.29 MB/s** | **164 MB/s** |
| xz -9 | 6.51x | 0.82 MB/s | 64 MB/s |
| zstd -19 | 5.20x | 1.07 MB/s | 119 MB/s |
| gzip -9 | 3.42x | 8.46 MB/s | 87 MB/s |
| lz4 -12 | 2.94x | 18.57 MB/s | 970 MB/s |

In other words: 29–43% below the strong compressors on ratio, but 4–5x faster
to compress, and the fastest to decompress among them. It is for "read one
block at random", not for "make the archive as small as possible".

## Building

```sh
make            # cpenc.o / bwt.o / cpool.o / libcpool_hook.so
make test       # 5 suites
make bench      # benchmark tools
```

Needs C11 and pthread. x86-64 builds with `-march=x86-64-v2` by default; the
AVX2 kernel needs `make ARCHFLAGS=-march=x86-64-v3` (see "Known limits").

## API

```c
#include "cpenc.h"

size_t cap = cpenc_max_compressed_size(n);      /* upper bound, all overhead included */
int  cpenc_compress(src, n, dst, &cap, order);  /* 0 on success; dst_len returns the real size */
int  cpenc_decompress(src, cap, dst, &n);       /* 0 on success; dst_len: capacity in, real size out */
```

- `order` (0–4): search effort. 0 is fastest, 4 gives the best ratio; libcfs
  uses 2.
- **Blocks are independent**: decompressing block k needs no earlier block.
  That is where random access comes from.
- **Deterministic**: the same input and the same `order` produce identical
  bytes (a test guards this).
- One block is capped at `CPENC_MAX_BLOCK` = 256MB; a decoder rejects anything
  larger.

## Format

```
[12-byte header: magic u16 | version u8 | stride u8 | orig u32 | comp u32][payload]
```

`stride` says which model this block used:

| stride | Name | Purpose |
|---|---|---|
| 0xFF | RAW | stored as-is (when compressing does not pay) |
| 0xFE | HEX | hex text → nibble folding |
| 0xFD | SEG | internal segmentation, each segment picks its own model |
| 0xFC | BWT | reorder, then delta/run coding |
| 0xFB / 0xFA | REF / REF2 | repeated strings as references (REF2 is a columnar token stream) |

Within a group the encoder has several tags (constant, LUT, Huffman, run
length, constant plus exceptions, bit planes); the details are in the comments
in `cpenc.c`. **The decoder trusts no length in the header**: out-of-range
values, wraparound, nesting depth and segment count all have bounds, and the
fuzzers exist to keep them honest.

## Testing

```sh
make test       # 5 suites: boundary sizes / special patterns / random round-trip / recompression / RLE
make fuzz SECS=300      # decoder: malformed streams, UBSan, guard pages
make fuzz-rt SECS=300   # encoder round-trip: what it produces must decode, byte for byte
```

The two fuzz targets ask **different questions**, and both are necessary:

- `fuzz`: can a hostile stream be decoded safely (bounds, undefined behaviour)?
- `fuzz-rt`: **does the encoder's own output decode?** That target has caught
  two real bugs — an unaligned header write and the hex layer dropping case —
  neither of which a decode-side target can see even in principle.

Needs `compiler-rt` (`apk add compiler-rt`). ASan does not run in this
project's development environment (PRoot), so the targets use UBSan plus
`PROT_NONE` guard pages.

## Known limits

- **The hex layer has a hard gate**: it engages only at ≥90% hex characters,
  <1% foreign characters, and uniform case. A dump with `0x` prefixes usually
  does **not** take this path (the `x` counts as a foreign character).
- **The AVX2 kernel barely shows end to end**: 3.44x on the kernel itself (1MB
  group, nz=16, 4 bits per symbol), but real data lands 100% on REF2/SEG
  blocks, where the time goes to reference copies rather than LUT expansion
  (94MB tar: 524 → 529 MB/s). It is for LUT-shaped blocks, not a general
  speedup.
- **Compression is slow**: 4–14 MB/s single-threaded. libcfs spreads that over
  8MB blocks, and the `cpool` layer adds a thread pool; do not expect
  single-threaded numbers to compete with gzip on speed.
- **No encryption, no checksums**: there is no checksum in a block. Damage
  detection relies on libcfs's `verify` (which re-decodes) or an external
  tool; the bytes themselves carry no integrity protection.
- The ratio knobs (`CPENC_SEG_SIZE`, the reference-layer thresholds) are
  empirical values from measurement, not derived.

## License

GNU General Public License v3.0 (`LICENSE`), with one **section 7(b)
additional term** (`NOTICE`): the copyright notice and the name "JXLAB" must
be preserved with the program, and a derived work must keep it alongside the
notices for its own changes.

The same license as [libcfs](../libcfs); the two link against each other, so
the terms have to match.
