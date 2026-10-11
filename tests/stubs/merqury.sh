#!/usr/bin/env bash
# Test stub for merqury.sh <reads.meryl> <asm.fa> <out_prefix>; writes
# outputs in the documented Merqury formats.
set -euo pipefail
asm=$(basename "$2"); asm=${asm%.fa}
printf '%s\t120\t1000000\t45.2\t3.0e-05\n' "$asm" > "$3.qv"
printf '%s\tall\t990000\t1000000\t99.0\n' "$asm" > "$3.completeness.stats"
