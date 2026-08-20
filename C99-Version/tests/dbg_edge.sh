#!/bin/bash
cd /mnt/d/Axioms/Axiom-Zspace-CODE/C99-Version
GEN=$(ls generated/purce_src_*.c | tr '\n' ' ')
gcc -std=c99 -O2 -fopenmp -Igenerated -Isrc -o /tmp/dbg_edge tests/dbg_edge.c \
    src/zspace_ingestion.c src/zspace_core.c $GEN -lm 2>&1 | grep -c 'error'
for w in 51 101 501 1001 1621; do
  echo "--- window $w ---"
  OMP_NUM_THREADS=1 /tmp/dbg_edge $w
done