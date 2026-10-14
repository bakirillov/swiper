#!/bin/bash

# This script runs Algorithm 1 and passive learning for Figure 2a, 2d and 2e

GS="../data/processed/GS_6nt.fasta"
ORIGINAL_DATA="../data/external/41586_2015_BFnature14592_MOESM46_ESM.xlsx"

QSIZE=$1
SETUP=$2
NROUNDS=$3
OUTPUT=$4


LIBRARIES=("library1" "library2")
SEEDS=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10")
SETUPS=("active" "random")
POOLS=("complete_pool" "golden_set")
ORTHOLOG_LIST=('wt_SpCas9' 'SpCas9_VQR' 'SpCas9_EQR' 'SpCas9_VRER' 'SpCas9_D1135E' 'St1Cas9_20nt' 'St1Cas9_21nt' 'SaCas9_21nt' 'SaCas9_23nt')



for ortholog in "${ORTHOLOG_LIST[@]}"; do
    for library in "${LIBRARIES[@]}"; do
        for pool in "${POOLS[@]}"; do
            for seed in "${SEEDS[@]}"; do

                RESULT="${OUTPUT}/${ortholog}/${library}/${pool}/${seed}/${SETUP}/${NROUNDS}"
                echo ""
                echo ${RESULT}
                echo ""
                
                if [ ${pool} = "complete_pool" ]; then
                    # If we are running on complete pool, then no pool file is needed
                    POOLFILE="None"
                else
                    # Use Golden Set pool file
                    POOLFILE=${GS}
                fi
                
                # Run Algorithm 1
                bash F2_Algorithm1.sh \
                    "${library}" \
                    "${pool}" \
                    "${ortholog}" \
                    "${QSIZE}" \
                    "${RESULT}" \
                    "${seed}" \
                    "${NROUNDS}" \
                    "${SETUP}" \
                    "${POOLFILE}"

                FIXED_LIB=$(echo "${library}" | sed 's/\([a-zA-Z]\)\([0-9]\)/\1 \2/') 
                echo "${FIXED_LIB}"

                # Run passive learning with the same setup as Algorithm 1
                python ../src/passive_learning.py -s ${seed} -p None -e ${RESULT}/${QSIZE}/result.xlsx -o ${RESULT}/${seed}.${QSIZE}.passive.tsv -d "${ORIGINAL_DATA};${FIXED_LIB};${ortholog}_PPDV"
            done
        done
    done
done