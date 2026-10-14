#!/bin/bash


LIBRARY=$1 # library1 or library2
POOL=$2 # complete_pool or golden_set
ORTHOLOG=$3 # wt Cas9 etc
QSIZE=$4 # query size (default 18)
OUTPUT=$5 # output path
SEED=$6 # random number seed
POOLFILE=$7 # file with pool data
SETUP=$8 # active learning or random

FILLING="../data/interim/fig2/fillings/${LIBRARY}/${POOL}/${ORTHOLOG}_PPDV.tsv"

if [ ${SETUP} = "active" ]; then
    python ../src/swiper_cli.py -d "${OUTPUT}/result.xlsx" -s ${SEED} -q ${QSIZE} -o "None" -p ${POOLFILE} -t ${OUTPUT}
else
    python ../src/swiper_cli.py -d "${OUTPUT}/result.xlsx" -s ${SEED} -q ${QSIZE} -o "None" -p ${POOLFILE} -r -t ${OUTPUT}
fi

# fill_from_data.py imitates the experiment and fills result.xlsx with known data from the experiment
python ../src/fill_from_data.py -d "${OUTPUT}/result.xlsx" -e ${FILLING} -o ${OUTPUT}