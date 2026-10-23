#!/bin/bash
# progress of the all-atom assembly runs: last line of obs.dat of each run
cd "$(dirname "$0")"
for d in runs/aa_assembly_*; do
	tail -1 $d/obs.dat | awk -v d=$(basename $d) '$1 ~ /^[0-9]/ {printf "%-24s step %9d  clusters %2d  largest %2d  parallel %2d (in register %2d)  antiparallel %2d (in register %2d)  S %.2f\n",d,$1,$2,$3,$6,$8,$7,$9,$12}'
	grep -q "^Final" $d/e.log && echo "    finished: $(grep ^Final $d/e.log)"
done
