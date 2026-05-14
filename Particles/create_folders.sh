#!/bin/bash

particles=("proton" "iron" "helium" "oxygen")

for p in "${particles[@]}"; do
    
    for e in $(seq 16.0 0.1 18.0); do
        energy_dir="lgE_${e}"
        
        for a in $(seq 0.0 0.1 0.9); do
            angle_dir="sin2_${a}"
            
            path="${p}/${energy_dir}/${angle_dir}"
            
            mkdir -p "$path"
        done
    done
done
