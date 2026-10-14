"""This file provides the code to run training for learning curves.
This file is to be ran as a standalone script.
"""

import os
import sys
import json
import argparse
import numpy as np
import pandas as pd
import os.path as op
from tqdm import tqdm
import swiper_cli as swiper
import matplotlib.pyplot as plt
from passive_learning import load_PPDV_data, get_Xy
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, balanced_accuracy_score, matthews_corrcoef, average_precision_score, roc_auc_score


def train_model(df, y_column, train_ix, test_ix, length=6):
    """
    Trains a Random Forest model on the provided DataFrame and calculates evaluation metrics.

    This function prepares the training and test sets, fits the Random Forest model, 
    and computes classification metrics.

    Parameters
    ----------
        df (pandas.DataFrame): The input DataFrame.
        y_column (str): The name of the target column.
        train_ix (np.array(int)): indices of training set.
        test_ix (np.array(int)): indices of test set.
        length (int, optional): Length of the sequences (e.g., 6). Defaults to 6.

    Returns
    -------
        tuple: A tuple containing:
            - y_hat (numpy.ndarray): Predicted labels for the test set.
            - metrics (list): A list of calculated metrics: accuracy, balanced accuracy, 
                and Matthews correlation coefficient.
    """
    X, y = get_Xy(df, y_column)
    rf = swiper.make_RF()
    rf.fit(X[train_ix], y[train_ix])
    y_hat = rf.predict(X[test_ix])
    y_hat_proba = rf.predict_proba(X[test_ix])
    metrics = [
        accuracy_score(y_pred=y_hat, y_true=y[test_ix]),
        balanced_accuracy_score(y_pred=y_hat, y_true=y[test_ix]),
        matthews_corrcoef(y_pred=y_hat, y_true=y[test_ix])
    ]
    return(y_hat, metrics)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d", "--dataset",
        dest="dataset",
        action="store", 
        help="set the dataset with NGS results", 
        default="initial.xlsx;library_1;wt_Cas9_PPDV",
    )
    parser.add_argument(
        "-t", "--tts",
        dest="tts",
        action="store", 
        help="set the training set size", 
        default="0.1",
    )
    parser.add_argument(
        "-s", "--seed",
        dest="seed",
        action="store", 
        help="set the seed for random number generation", 
        default="10",
    )
    parser.add_argument(
        "-o", "--output",
        dest="output",
        action="store", 
        help="set the output data filename", 
        default="passive.tsv",
    )
    args = parser.parse_args()
    seed = int(args.seed)
    np.random.choice(seed)
    PATH, sheet_name, ortholog = args.dataset.split(";")
    sheet_name = sheet_name.replace("_", " ")
    ortholog = ortholog.replace("_", " ")
    data = load_PPDV_data(excel_file=PATH, sheet_name=sheet_name)
    print(data.shape)
    traintest_ix, test_ix = train_test_split(np.arange(data.shape[0]), test_size=0.25, random_state=seed)
    train_ix, _ = train_test_split(traintest_ix, train_size=float(args.tts), random_state=seed)
    print(train_ix.shape, test_ix.shape)
    _, metrics = train_model(data, ortholog, train_ix, test_ix, length=6)
    os.makedirs(op.split(args.output)[0], exist_ok=True)
    with open(args.output.replace(".tsv", ".metrics.tsv"), "w") as oh:
        oh.write("\t".join([str(a) for a in metrics])+"\n")