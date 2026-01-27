#!/bin/bash

primary_particle="proton"
logE="18.0"
sin2theta="0.4"

# Construct paths
SRC="rawdata/${primary_particle}/lgE_${logE}/sin2_${sin2theta}"
DEST="DATfiles/${primary_particle}/lgE_${logE}/sin2_${sin2theta}"

# Create destination directory (including parent directories if not exist)
mkdir -p "$DEST"

# Loop over run-number directories inside SRC
for dir in "$SRC"/[0-9][0-9][0-9][0-9][0-9][0-9]; do
    [ -d "$dir" ] || continue  # skip if not a directory

    run=$(basename "$dir")
    file="$dir/DAT$run"

    if [ -f "$file" ]; then
        cp "$file" "$DEST/"
    else
        echo "Warning: $file not found"
    fi
done
