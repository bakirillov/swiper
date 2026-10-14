#!/bin/bash

# This script computes the query size selection curve for Figure 2
# Takes LONG time to compute, you are better off running the scripts in parallel

for ((i=3; i<=35; i++)); do
    bash ./F2ade.sh ${i} active 10 ../data/interim/fig2/
done

for num in 40 50 100; do
    bash ./F2ade.sh ${num} active 10 ../data/interim/fig2/
done

