#!/bin/bash

# This script computes the learning curve for passive learning (panels b and c of Figure 2)

OUTPUT=$1

ORIGINAL_DATA="../data/external/41586_2015_BFnature14592_MOESM46_ESM.xlsx"
TRAININGSETSIZES=("0.01" "0.05" "0.1" "0.2" "0.3" "0.4" "0.5" "0.6" "0.7" "0.8" "0.9" "0.95" "0.99" "1.0")
LIBRARIES=("library1" "library2")
SEEDS=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10")
ORTHOLOG_LIST=('wt_SpCas9' 'SpCas9_VQR' 'SpCas9_EQR' 'SpCas9_VRER' 'SpCas9_D1135E' 'St1Cas9_20nt' 'St1Cas9_21nt' 'SaCas9_21nt' 'SaCas9_23nt')


for ortholog in "${ORTHOLOG_LIST[@]}"; do
    for library in "${LIBRARIES[@]}"; do
        for seed in "${SEEDS[@]}"; do
            for tts in "${TRAININGSETSIZES[@]}"; do
                RESULT="${OUTPUT}/${ortholog}/${library}/passive/${tts}_${seed}"
                echo ""
                echo ${RESULT}
                echo ""

                POOLFILE="/tmp/passive_${ortholog}_${library}_${seed}_${tts}.fasta"
                FIXED_LIB=$(echo "${library}" | sed 's/\([a-zA-Z]\)\([0-9]\)/\1 \2/') 
                echo "${FIXED_LIB}"
                
                python ../src/learning_curve.py -s ${seed} --tts ${tts} -o ${RESULT}.tsv -d "${ORIGINAL_DATA};${FIXED_LIB};${ortholog}_PPDV"
            done
        done
    done
done