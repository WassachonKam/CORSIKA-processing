#!/bin/bash

particles=("proton" "iron" "helium" "oxygen")

for p in "${particles[@]}"; do
    
    for e in $(seq 16.0 0.1 18.0); do
        energy_dir="lgE_${e}"
        path="${p}/${energy_dir}"
        
        mkdir -p "$path"

    done
done
