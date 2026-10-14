#!/bin/bash

# This script gathers the results of data for Figure 2
# Takes LONG time to compute, you are better off running the scripts in parallel

for ((i=3; i<=35; i++)); do
    python ../src/gather_results.py -q ${i} -o ../data/processed/F2/${i}.tsv
done

for num in 40 50 100; do
    python ../src/gather_results.py -q ${num} -o ../data/processed/F2/${num}.tsv
done