"""This file provides the code to gather training results for later processing.
This file is to be ran as a standalone script.
"""


import os
import re
import sys
import json
import einops
import argparse
import itertools
import numpy as np
import pandas as pd
import pickle as pkl
import os.path as op
from tqdm import tqdm
from Bio.Seq import Seq
import Levenshtein as lev
from Bio import pairwise2
import multiprocessing as mp
from itertools import product
from Bio import SeqIO, motifs
from itertools import product
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.colors import LinearSegmentedColormap
from sklearn.metrics import matthews_corrcoef, accuracy_score, balanced_accuracy_score

import warnings
warnings.filterwarnings('ignore')

try:
    from . import swiper_cli as swiper
except ImportError:
    import swiper_cli as swiper
    
try:
    from . import convergence_test as convergence
except ImportError:
    import convergence_test as convergence
    
    
random_query_strategy = swiper.random_query_strategy
    

class Result():
    """
    Container for managing CRISPR PAM screening result paths and data retrieval.

    Stores experiment parameters, file paths, and consensus sequences.
    Provides methods to extract passive and active learning metrics.

    Args:
        n_queries (str): Number of queries.
        n_round (str): Number of rounds.
        base (str): Base directory path.

    Methods:
        make_path(ortholog, library, pool, seed, setup):
            Constructs file paths for test, metrics, and model files.
        passive():
            Reads passive metrics from file.
        active():
            Loads model, predicts, and calculates metrics.
        line():
            Assembles the result row for the output DataFrame.
    """

    
    def __init__(
        self, n_queries, n_round, base="../data/interim/fig2/ORTHOLOG/LIBRARY/POOL/SEED/SETUP/10/"
    ):
        self.BASE = base
        self.TEST = "SEED.NQUERIES.passive.test.tsv"
        self.PASSIVE_METRICS = "SEED.NQUERIES.passive.metrics.tsv"
        self.MODEL = "NQUERIES/model_NROUND.pkl"
        self.NQUERIES = str(n_queries)
        self.NROUND = str(n_round)
        self.columns = [
            "ortholog", "library", "pool", "seed", "setup", "n_queries", "n_round",
            "passive_acc", "passive_balaccs", "passive_mcc", 
            "active_acc", "active_balaccs", "active_mcc", "pam", "score"
        ]
        self.consensi = {
            'library1': {'wt_SpCas9': 'SGGGMN',
                'SpCas9_VQR': 'NGNGSN',
                'SpCas9_EQR': 'SGAKGM',
                'SpCas9_VRER': 'DGCGGC',
                'SpCas9_D1135E': 'SGGGNN',
                'St1Cas9_20nt': 'NNANAA',
                'St1Cas9_21nt': 'SNRGAA',
                'SaCas9_21nt': 'SMGART',
                'SaCas9_23nt': 'SMGRRW'},
            'library2': {'wt_SpCas9': 'SGGGNN',
                'SpCas9_VQR': 'NGNGNN',
                'SpCas9_EQR': 'NGAKGN',
                'SpCas9_VRER': 'DGCGGC',
                'SpCas9_D1135E': 'SGGGNN',
                'St1Cas9_20nt': 'CWARAA',
                'St1Cas9_21nt': 'MNRGAA',
                'SaCas9_21nt': 'SNGRNN',
                'SaCas9_23nt': 'SNGRNN'}
        }

    def make_path(self, ortholog, library, pool, seed, setup):
        """
        Constructs file paths by substituting placeholders in the base directory
        with experiment parameters.

        This method updates the instance attributes (`self.BASE`, `self.TEST`, etc.)
        to reflect the constructed paths for test, passive metrics, and model files.

        Parameters
        ----------
            ortholog (str): The identifier for the ortholog.
            library (str): The identifier for the library.
            pool (str): The identifier for the pool.
            seed (str): The identifier for the seed.
            setup (str): The identifier for the setup.

        Returns
        -------
            None
        """
        self.BASE = self.BASE.replace(
            "ORTHOLOG", ortholog
        ).replace("LIBRARY", library).replace("POOL", pool).replace("SEED", seed).replace("SETUP", setup)
        self.TEST = self.BASE+self.TEST.replace("SEED", seed).replace("NQUERIES", self.NQUERIES)
        self.PASSIVE_METRICS = self.BASE+self.PASSIVE_METRICS.replace(
            "SEED", seed
        ).replace("NQUERIES", self.NQUERIES)
        self.MODEL = self.BASE+self.MODEL.replace("NROUND", self.NROUND).replace("NQUERIES", self.NQUERIES)
        self.setups = [ortholog, library, pool, seed, setup, self.NQUERIES, self.NROUND]

    def passive(self):
        """
        Reads passive metrics from a CSV file.

        Attempts to load a tab-separated CSV file located at `self.PASSIVE_METRICS`.
        If the file is missing or cannot be read, it returns a placeholder DataFrame
        with NaN values for columns 0, 1, and 2.

        Returns:
            list: A list of lists representing the rows of the data.
        """
        try:
            data = pd.read_csv(self.PASSIVE_METRICS, sep="\t", header=None, index_col=None)
            data.columns = [0,1,2]
        except:
            data = pd.DataFrame(
                {0: [np.nan], 1: [np.nan], 2: [np.nan]}
            )
        finally:
            return(data.values.tolist())
    
    def active(self):
        """
        Loads the trained model and test data to evaluate predictions and calculate consensus scores.

        This method loads the model pickle and test data, preprocesses the PAM features,
        performs predictions, and calculates classification metrics (accuracy, balanced accuracy, MCC).
        It also identifies the consensus sequence for the active selection and computes its score.

        Returns:
            list: A list containing:
                - Model evaluation metrics for predictions (Accuracy, Balanced Accuracy, MCC).
                - The active consensus sequence.
                - The consensus score.
        """
        test = pd.read_csv(self.TEST, sep="\t", index_col=0)
        current_X = np.stack(test["PAM"].apply(swiper.onehot).values).astype(np.float32)
        current_Y = test["Outcome"].values        
        with open(self.MODEL, "rb") as ih:
            current_model = pkl.load(ih)
            y_hat = current_model.predict(current_X)
            classified, _ = swiper.classify_all(current_model, 6)
            try:
                current_pam_attempt = convergence.get_consensus(
                    classified[classified["Y_class"] == "+"]["PAM"].values
                )
            except ValueError:
                current_pam_attempt = "N"*6
            finally:
                current_pam_score = convergence.consensus_score(
                    self.consensi[self.setups[1]][self.setups[0]], current_pam_attempt
                )
            active = [
                accuracy_score(y_pred=y_hat, y_true=current_Y),
                balanced_accuracy_score(y_pred=y_hat, y_true=current_Y),
                matthews_corrcoef(y_pred=y_hat, y_true=current_Y),
                current_pam_attempt,
                current_pam_score
            ]
        return(active)

    def line(self):
        """
        Assembles the final row structure for the output DataFrame by combining configuration parameters,
        passive metrics, and active prediction results.

        Prepares a flat list starting with setup parameters, extends it with passive metrics, 
        and then extends it with active prediction metrics. Converts the resulting list into a 
        pandas DataFrame with the predefined column names.

        Returns:
            pandas.DataFrame: The assembled result row.
        """
        line = self.setups[:]
        line.extend(*self.passive())
        line.extend(self.active())
        out = pd.DataFrame(np.array(line)).T
        out.columns = self.columns
        return(out)



def process_combination(combination):
    """
    Processes a single combination of (directory, b, query).
    Returns the line (pandas Series/DataFrame) or None on error.
    """
    a,b,c = combination
    r = Result(c, str(b))
    r.make_path(*a)
    line = r.line()
    return line



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-q", "--qsize",
        dest="qsize",
        action="store", 
        help="set the query size", 
        default="18",
    )
    parser.add_argument(
        "-o", "--output",
        dest="output",
        action="store", 
        help="set the output data filename", 
        default="18.tsv",
    )
    args = parser.parse_args()
    LIBRARIES = ["library1", "library2"]
    SEEDS = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10"]
    ORTHOLOG_LIST = [
        'wt_SpCas9', 'SpCas9_VQR', 'SpCas9_EQR', 'SpCas9_VRER', 'SpCas9_D1135E', 
        'St1Cas9_20nt', 'St1Cas9_21nt', 'SaCas9_21nt', 'SaCas9_23nt'
    ]
    SETUPS = ["active", "random"]
    SWITCHES = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10"]
    NQUERIES_ = [str(a) for a in np.arange(1,31)]+["35", "40", "50", "100"]
    POOLS = ["complete_pool", "golden_set"]
    QUERY_SIZE = str(args.qsize)
    print(QUERY_SIZE)
    fn = str(args.output)
#    fn = "../data/processed/F2/"+QUERY_SIZE+".tsv"
    result = []
    DIRECTORIES = itertools.product(ORTHOLOG_LIST, LIBRARIES, POOLS+SWITCHES, SEEDS, SETUPS)
    combinations = list(product(DIRECTORIES, np.arange(1, 11), [str(QUERY_SIZE)]))
    with open(fn+".failed", "a") as oh:
        for a in tqdm(combinations):
            line = process_combination(a)
   #         print(a, line)
            if line.shape[0] > 0:
                result.append(line)
            else:
                oh.write(str(line)+";"+str(a)+"\n")
    result = pd.concat(result)
    result.to_csv(fn, sep="\t")
    print(fn+" DONE")