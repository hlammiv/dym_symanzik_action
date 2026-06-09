#!/bin/bash
# Level-2 parallel 2D-LLR: one sequential seed pass (drive-in + dump configs),
# then NC parallel solve workers over disjoint cell ranges, then combine.
# Reproduces `llr2d` full-mode output but parallelizes the (expensive) RM over
# cells -- use it for larger lattices where a single full run is slow.
#
# Usage:
#   llr2d_parallel.sh <outdir> <NC> <group> <D> <Nt> <Nx> \
#       <E1top> <E1bot> <step1> <hw1> <c0> <c1> <c2> <step2> <hw2> <nperp> \
#       <a0> <seed> <K> <NRM>
# Produces <outdir>/combined.out  (feed to scripts/llr2d_reconstruct.py).
set -e
OUT=$1; NC=$2; shift 2
GRP=$1; D=$2; Nt=$3; Nx=$4; shift 4
E1top=$1; E1bot=$2; step1=$3; hw1=$4; c0=$5; c1=$6; c2=$7
step2=$8; hw2=$9; nperp=${10}; a0=${11}; seed=${12}; K=${13}; NRM=${14}
LLR2D=${LLR2D:-./llr2d}
mkdir -p "$OUT/cfg"

echo "[seed] sequential drive-in -> $OUT/cfg"
"$LLR2D" seed "$GRP" "$D" "$Nt" "$Nx" "$E1top" "$E1bot" "$step1" "$hw1" \
    "$c0" "$c1" "$c2" "$step2" "$hw2" "$nperp" "$a0" "$seed" "$OUT/cfg" \
    > "$OUT/seed.log" 2> "$OUT/seed.err"
grep '^CELL:' "$OUT/seed.log" > "$OUT/manifest.txt"
N=$(wc -l < "$OUT/manifest.txt")
echo "[seed] $N cells"

echo "[solve] $NC workers"
per=$(( (N + NC - 1) / NC ))
for w in $(seq 0 $((NC-1))); do
  lo=$((w*per)); hi=$(((w+1)*per))
  [ "$lo" -ge "$N" ] && break
  OMP_NUM_THREADS=1 "$LLR2D" solve "$GRP" "$D" "$Nt" "$Nx" "$OUT/manifest.txt" \
      "$lo" "$hi" "$a0" "$K" "$NRM" "$((seed+1000+w))" \
      > "$OUT/solve_$w.out" 2>> "$OUT/solve.err" &
done
wait

grep '^LLR2D' "$OUT/seed.log" > "$OUT/combined.out"
cat "$OUT"/solve_*.out | grep '^ANE2:' >> "$OUT/combined.out"
echo "[done] $(grep -c '^ANE2:' "$OUT/combined.out") cells -> $OUT/combined.out"
