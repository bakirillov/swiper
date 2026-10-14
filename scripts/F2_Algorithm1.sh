#!/bin/bash

# This script implements Algorithm 1 from the paper


LIBRARY=$1 # library1 or library2
POOL=$2 # complete_pool or golden_set
ORTHOLOG=$3 # wt Cas9 etc
QSIZE=$4 # query size (default 18)
OUTPUT=$5 # output path
SEED=$6 # random number seed
ROUNDS=$7 # number of rounds
SETUP=$8 # active or random
POOLFILE=$9 # file with pool data, e.g. ../data/processed/GS_6nt.fasta, None if entire pool is used

FILLING="../data/interim/fig2/fillings/${LIBRARY}/${POOL}/${ORTHOLOG}_PPDV.tsv"

echo ${FILLING}

RESULT="${OUTPUT}/${QSIZE}"


# Create the result.xlsx to store the experimental state
# The scripts save state to result.xlsx and rewrite it
python ../src/initialize.py -d ${FILLING} -q ${QSIZE} -o ${RESULT} -s ${SEED}


for ((i=1; i<=ROUNDS; i++)); do
    # Run single round of experiment
    echo "Round $i"
    bash F2_Single_Round.sh ${LIBRARY} ${POOL} ${ORTHOLOG} ${QSIZE} ${RESULT} ${SEED} ${POOLFILE} ${SETUP}
done


