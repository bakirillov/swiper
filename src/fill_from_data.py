"""This file provides the code to fill the input data for Proof-of-Concept study.
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
        default="initial.xlsx",
    )    
    parser.add_argument(
        "-e", "--experiments",
        dest="experiments",
        action="store", 
        help="set the experimental data to use to fill the data with the experiments"
    )
    parser.add_argument(
        "-o", "--output",
        dest="output",
        action="store", 
        help="set the output data directory"
    )
    args = parser.parse_args()
    xls = pd.ExcelFile(args.dataset)
    history = pd.read_excel(xls, None)
    history = {a: history[a][["#", "PAM", "Outcome"]] for a in history}
    experiments = pd.read_csv(args.experiments, sep="\t", index_col=0)
    print(np.unique(experiments["Outcome"].values, return_counts = True))
    latest = str(list(sorted([int(a) for a in list(history.keys())]))[-1])
    latest_df = history[latest]
    b = np.isin(experiments["PAM"], latest_df["PAM"])
    small = experiments.loc[b]
    latest_df["Outcome"] = [str(small[small.PAM == a].Outcome.values[0]) for a in latest_df["PAM"]]
    history[latest] = latest_df
    swiper.make_history(history, args.output)