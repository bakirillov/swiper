"""This file provides the code to initialize an experiment.
This file is to be ran as a standalone script.
"""


import os
import json
import shutil
import argparse
import itertools
import numpy as np
import pandas as pd
import pickle as pkl
import os.path as op
from tqdm import tqdm
import swiper_cli as swiper


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d", "--dataset",
        dest="dataset",
        action="store", 
        help="set the dataset to fill with the experiments", 
    )
    parser.add_argument(
        "-q", "--query-size",
        dest="qsize",
        action="store", 
        help="18", 
    )
    parser.add_argument(
        "-o", "--output",
        dest="output",
        action="store", 
        help="set the output data directory"
    )
    parser.add_argument(
        "-s", "--seed",
        dest="seed",
        action="store", 
        help="set the seed for random number generation", 
        default="10",
    )
    args = parser.parse_args()
    outdir = args.output
    qsize = int(args.qsize)
    seed = int(args.seed)
    np.random.seed(seed)
    library = pd.read_csv(args.dataset, sep="\t")
    sheet = library.sample(qsize, random_state=seed)
    sheet.columns = ["#", "PAM", "Outcome"]
    #print(current_sample)
    #make_new_sheet(Ns, PAMs)
    os.makedirs(outdir, exist_ok=True)
    swiper.make_history({"1": sheet}, outdir)