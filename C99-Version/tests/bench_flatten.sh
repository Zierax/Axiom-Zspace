#!/bin/bash
cd /mnt/d/Axioms/Axiom-Zspace-CODE/C99-Version
GEN=$(ls generated/purce_src_*.c | tr '\n' ' ')
gcc -std=c99 -O2 -fopenmp -Igenerated -Isrc -o /tmp/bench_flatten tests/bench_flatten.c \
    src/zspace_ingestion.c src/zspace_core.c $GEN -lm 2>&1 | grep -c 'error'
echo '--- 1 thr ---'
OMP_NUM_THREADS=1 /tmp/bench_flatten 91566 1567
echo '--- 16 thr ---'
OMP_NUM_THREADS=16 /tmp/bench_flatten 91566 1567
echo '--- 1 thr small window ---'
OMP_NUM_THREADS=1 /tmp/bench_flatten 91566 51
echo '--- 16 thr small window ---'
OMP_NUM_THREADS=16 /tmp/bench_flatten 91566 51