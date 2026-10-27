#!/usr/bin/env bash
# Generate a parameter list: one line per (primary, lgE, third column) combination.
# Usage: ./make_params.sh [output_file]   (default: params_all.txt)

OUT="${1:-params_all.txt}"

PRIMARIES=(proton iron helium oxygen)

# lgE from 16.0 to 18.0 in steps of 0.1 (stored as integers x10 to avoid float drift)
LGE_MIN=160
LGE_MAX=180

# Third column from 0.0 to 0.9 in steps of 0.1 (stored as integers x10)
VAL_MIN=2
VAL_MAX=3

: > "$OUT"
for p in "${PRIMARIES[@]}"; do
    for ((e = LGE_MIN; e <= LGE_MAX; e++)); do
        for ((v = VAL_MIN; v <= VAL_MAX; v++)); do
            printf "%s lgE_%d.%d %d.%d\n" "$p" $((e / 10)) $((e % 10)) $((v / 10)) $((v % 10)) >> "$OUT"
        done
    done
done

echo "Wrote $(wc -l < "$OUT") lines to $OUT"