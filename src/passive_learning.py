"""This file provides the code to run training for passive learning baseline.
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
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, balanced_accuracy_score, matthews_corrcoef, average_precision_score, roc_auc_score

try:
    from . import swiper_cli as swiper
except ImportError:
    import swiper_cli as swiper



def load_PPDV_data(excel_file, sheet_name="library 1"):
    """
    Loads PPDV data from an Excel file.
    
    Reads the specified sheet, drops the header row, renames columns to standard names,
    and calculates relative abundances (PPDV) for various Cas9 variants by dividing
    their counts by the 'library 1 frequency'.
    
    Pre-Post Depletion Value is defined as the ratio of read numbers before and
    after cleavage, normalized on library size. For all experiments, following the original
    paper [29], we used PPDV = 0.2 as a threshold to binarize the labels 
    into successful (PPDV ≤ 0.2) and unsuccessful (PPDV > 0.2) cleavage. 

    Parameters
    ----------
    excel_file : str
        Path to the PPDV file.
    sheet_name : str
        Name of the shit to load.

    Returns
    -------
    pd.DataFrame
        A dataframe that contains PPDV data.
    
    """
    df = pd.read_excel(excel_file, engine="openpyxl", sheet_name=sheet_name)
    cols = df.iloc[0]
    df = df.iloc[1:]
    new_cols = [
        "6nt PAM", 
        "library 1 count", "library 1 frequency",
        "dead SpCas9 count", "dead SpCas9 frequency",
        "wt SpCas9 count", "wt SpCas9 frequency",
        "SpCas9 VQR count", "SpCas9 VQR frequency",
        "SpCas9 EQR count", "SpCas9 EQR frequency",
        "SpCas9 VRER count", "SpCas9 VRER frequency",
        "SpCas9 D1135E count", "SpCas9 D1135E frequency",
        "St1Cas9 20nt count", "St1Cas9 20nt frequency",
        "St1Cas9 21nt count", "St1Cas9 21nt frequency",
        "SaCas9 21nt count", "SaCas9 21nt frequency",
        "SaCas9 23nt count", "SaCas9 23nt frequency",
    ]
    df.columns = new_cols
    freq_cols = list(filter(lambda x: "frequency" in x, df.columns))[1:]
    for a in freq_cols:
        newcol = a.replace("frequency", "PPDV")
        df[newcol] = df[a]/df['library 1 frequency']
    return(df)


def get_Xy(df, y_column, threshold=0.2):
    """
    Prepares feature (X) and target (y) arrays for a binary classification model.

    This function extracts the '6nt PAM' column (features) and the target column.
    It applies a threshold to binarize the target (1 if <= threshold, else 0) and
    one-hot encodes the features using the swiper.onehot function.
    
    Parameters
    ----------
        df (pandas.DataFrame): The input DataFrame.
        y_column (str): The name of the column to be used as the target variable.
        threshold (float, optional): The threshold value for classification.
            Defaults to 0.2.

    Returns
    -------
        tuple: A tuple containing:
            - X (numpy.ndarray): The feature matrix (one-hot encoded).
            - y (numpy.ndarray): The target labels (0 or 1).
    """
    Xdf = df[["6nt PAM", y_column]]
    y = Xdf[y_column].apply(lambda x: 1 if x <= threshold else 0).values
    X = np.stack(Xdf["6nt PAM"].apply(swiper.onehot).values)
    return(X,y)


def train_model(
    df, y_column, seed=1, length=6, pool_fn="None", train_size=100, threshold=0.2, stratify=False
):
    """
    Trains a Random Forest model on the provided DataFrame and calculates evaluation metrics.

    This function prepares the training and test sets, optionally filters the training data based on a
    sequence pool file, fits the Random Forest model, and computes classification metrics.

    Parameters
    ----------
        df (pandas.DataFrame): The input DataFrame.
        y_column (str): The name of the target column.
        seed (int, optional): Random seed for reproducibility. Defaults to 1.
        length (int, optional): Length of the sequences (e.g., 6). Defaults to 6.
        pool_fn (str, optional): Path to the sequence pool file. If not "None", the training set
            is filtered by this pool. Defaults to "None".
        train_size (int, optional): Number of samples to include in the training set. Defaults to 100.
        threshold (float, optional): Threshold value for binarizing the target. Defaults to 0.2.
        stratify (bool, optional): Whether to stratify the train-test split. Defaults to False.

    Returns
    -------
        tuple: A tuple containing:
            - y_hat (numpy.ndarray): Predicted labels for the test set.
            - metrics (list): A list of calculated metrics: accuracy, balanced accuracy, 
                and Matthews correlation coefficient.
            - test_ix (numpy.ndarray): Indices of the test set.
    """
    X, y = get_Xy(df, y_column)
    if stratify:
        train_ix, test_ix = train_test_split(np.arange(X.shape[0]), stratify=y, random_state=seed)
    else:
        train_ix, test_ix = train_test_split(np.arange(X.shape[0]), random_state=seed)
    if pool_fn != "None":
        pool = []
        with open(args.pool, "r") as ih:
            for a in ih.readlines():
                if a[0] != ">":
                    pool.append(a.replace("\n", ""))
        pool = np.array(pool)
        b_pool = np.isin(data.iloc[train_ix]["6nt PAM"], pool)
        train_ix = train_ix[b_pool]
        print(train_ix.shape)
    train_ix = train_ix[0:train_size]
    rf = swiper.make_RF()
    rf.fit(X[train_ix], y[train_ix])
    y_hat = rf.predict(X[test_ix])
    y_hat_proba = rf.predict_proba(X[test_ix])
    metrics = [
        accuracy_score(y_pred=y_hat, y_true=y[test_ix]),
        balanced_accuracy_score(y_pred=y_hat, y_true=y[test_ix]),
        matthews_corrcoef(y_pred=y_hat, y_true=y[test_ix])
    ]
    return(y_hat, metrics, test_ix)


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
        "-e", "--exclude",
        dest="exclude",
        action="store", 
        help="exclude the processed PAM candidates from the training pool", 
        default="exclude.xlsx",
    )
    parser.add_argument(
        "-p", "--pool",
        dest="pool",
        action="store", 
        help="set the pool to select PAM candidates from", 
        default="None",
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
    xls = pd.ExcelFile(args.exclude)
    history = pd.read_excel(xls, None)
    history = {a: history[a][["#", "PAM", "Outcome"]] for a in history}
    exclude = list(set(sum([history[a]["PAM"].values.tolist() for a in history], [])))
    print("Excluding ", len(exclude))
    b_ex = np.logical_not(np.isin(data["6nt PAM"], exclude))
    print(data.shape)
    data = data.loc[b_ex]
    print(data.shape)
    _, metrics, test_ix = train_model(data, ortholog, seed=seed, pool_fn=args.pool, train_size=int(data.shape[0]*0.75), stratify=False)
    with open(args.output.replace(".tsv", ".metrics.tsv"), "w") as oh:
        oh.write("\t".join([str(a) for a in metrics])+"\n")
    test_data = data.iloc[test_ix]
    _, y = get_Xy(test_data, ortholog, 0.2)
    pd.DataFrame(
        {
            "PAM": test_data["6nt PAM"].values, 
            "Outcome": ["-" if a == 0 else "+" for a in y]
        }
    ).to_csv(args.output.replace(".tsv", ".test.tsv"), sep="\t")
    
    
# REFERENCES:
# 29. Kleinstiver, B. P. et al. Engineered crispr-cas9 nucleases with altered pam specificities. Nature 523, 481–485 (2015).