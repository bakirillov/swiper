#!/bin/bash

# This script implements Algorithm 2 from the paper

LIBRARY=$1 # library1 or library2
SWITCH=$2 # round to switch on
ORTHOLOG=$3 # wt Cas9 etc
QSIZE=$4 # query size (default 18)
OUTPUT=$5 # output path
SEED=$6 # random number seed
ROUNDS=$7 # number of rounds
SETUP=$8 # active or random
POOLFILE=$9 # file with pool data, e.g. ../data/processed/GS_6nt.fasta, None if entire pool is used

FILLING="../data/interim/fig2/fillings/${LIBRARY}/golden_set/${ORTHOLOG}_PPDV.tsv"

echo ${FILLING}

RESULT="${OUTPUT}/${QSIZE}"

# Create the result.xlsx to store the experimental state
# The scripts save state to result.xlsx and rewrite it
python ../src/initialize.py -d ${FILLING} -q ${QSIZE} -o ${RESULT} -s ${SEED}


for ((i=1; i<=ROUNDS; i++)); do
    echo "Round $i"
    if (( SWITCH > i )); then
        # if switch round is not reached, run single round on the Golden Set pool
        bash F2_Single_Round.sh ${LIBRARY} golden_set ${ORTHOLOG} ${QSIZE} ${RESULT} ${SEED} ${POOLFILE} ${SETUP}
    else
        # if switch round is reached, run single round on the entire pool of all possible k-mers
        echo "SWITCHED"
        bash F2_Single_Round.sh ${LIBRARY} complete_pool ${ORTHOLOG} ${QSIZE} ${RESULT} ${SEED} None ${SETUP}
    fi
done


