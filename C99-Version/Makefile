# ═══════════════════════════════════════════════════════════════════════════
# C99-Version Makefile  (WSL gcc, -std=c99, no external deps)
# ═══════════════════════════════════════════════════════════════════════════
CC      := gcc
CFLAGS  := -std=c99 -O3 -march=native -mtune=native -flto -ffast-math -fno-math-errno -funroll-loops -DNDEBUG -Wall -Wextra -Wno-unused-parameter -fopenmp -fopenmp-simd -Igenerated -Isrc
GEN     := $(wildcard generated/purce_src_*.c)

all: bin/verify_kernels bin/zspace_card

bin/verify_kernels: tests/verify_kernels.c $(GEN) | bin
	$(CC) $(CFLAGS) -o $@ tests/verify_kernels.c $(GEN) -lm

bin/zspace_card: src/zspace_card.c src/zspace_core.c src/zspace_bls.c \
                 src/zspace_eph.c src/zspace_audit.c src/zspace_ingestion.c \
                 src/zspace_batch.c $(GEN) | bin
	$(CC) $(CFLAGS) -o $@ src/zspace_card.c src/zspace_core.c \
	                 src/zspace_bls.c src/zspace_eph.c src/zspace_audit.c \
	                 src/zspace_ingestion.c src/zspace_batch.c \
	                 $(GEN) -lm

bin:
	mkdir -p bin

clean:
	rm -rf bin

.PHONY: all clean