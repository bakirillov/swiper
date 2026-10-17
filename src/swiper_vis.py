"""
Script to visualize swiper results.
Takes a left.zip file (or pair of left.zip and right.zip) and generates an output.pdf file.
"""

import os
import sys
import argparse
import numpy as np
import os.path as op
import pickle as pkl
import matplotlib.pyplot as plt

try:
    from . import swiper_cli as swiper
except ImportError:
    import swiper_cli as swiper
    
try:
    from . import convergence_test as convergence
except ImportError:
    import convergence_test as convergence


HEADER_TEMPLATE = """

# EXPERIMENT_NAME

"""
    
ITER_TEMPLATE = """

## ITERNUMBER-th iteration

Predicted PAM consensus: consensus_string      

![Predicted PAM logo at ITERNUMBER-th iteration](consensus_logo)       

"""

CONVERGENCE_TEMPLATE = """
## Statistical test for convergence of PAM search

![CONVERGENCE](convergence_violin)     
"""
    


class Output():
    """
    A class to generate formatted output strings from a template file.

    Parameters
    ----------
    name (str): Name of the experiment.
    consensi (list): List containing consensus strings and logos.
    right (list, optional): List containing right consensus string and logo. Defaults to None.
    convergence (list, optional): List containing convergence details. Defaults to None.
    """
    
    def __init__(self, name, consensi, convergence=None):
        """
        Initializes the Output instance.

        Reads the template file, parses it, and performs string replacements
        using the provided arguments.

        Parameters
        ----------
        name (str): Name of the experiment.
        consensi (list): List containing consensus strings and logos.
        right (list, optional): List containing right consensus string and logo. Defaults to None.
        convergence (list, optional): List containing convergence details. Defaults to None.
        """
        self.OUTPUT = [
            HEADER_TEMPLATE.replace("EXPERIMENT_NAME", name)
        ]
        for i,(consensus, logo) in enumerate(consensi):
            self.OUTPUT.append(
                ITER_TEMPLATE.replace("consensus_string", consensus).replace(
                    "ITERNUMBER", str(i+1)
                ).replace(
                    "consensus_logo", logo
                )
            )
        if convergence:
            self.OUTPUT.append(
                CONVERGENCE_TEMPLATE.replace("convergence_violin", convergence[1])
            )
        
    def out(self):
        """
        Returns the formatted output string.

        Returns:
            str: The formatted string.
        """
        return("\n\n".join(self.OUTPUT))


class ModelResult():
    
    def __init__(self, pkl_path, n=7, n_value=0.1):
        """
        Class for visualization of results for a single test and a pair of tests.

        Parameters
        ----------
        pkl_path (str): path to model_N.pkl.
        n (int, optional): PAM candidate length used during classification
                           (passed directly to swiper.classify_all). Defaults to 7.
        n_value (float, optional): Scaling/weighting factor for the convergence score
                                   calculation. Defaults to 0.1.
        """
        self.n = n
        self.n_value = n_value
        self.path = pkl_path    

    @staticmethod
    def draw_logo(df, fn):
        """
        Generates a sequence logo image from the PAM column of a DataFrame.

        Parameters
        ----------
        df (pandas.DataFrame): Must contain a column named "PAM".
        fn (str): Absolute or relative filepath where the .png logo will be saved.
        """
        rnd = "".join(np.random.choice(np.arange(10), size=10).astype(str))
        fn_ = "/tmp/swiper_"+rnd+".fasta"
        swiper.save_fasta(fn_, df["PAM"])
        swiper.make_logo(fn_, fn)
    
    def load_model(self):
        """
        Loads a pickled classification model from the disk.

        Each protein and each active learning round has its own trained model.
        These models are the core classifiers used to predict which PAMs are cleaved.
        
        Returns
        -------
        object: Deserialized model.
        """
        with open(self.path, "rb") as ih:
            model = pkl.load(ih)
        return(model)

    def screen(self):
        """
        Runs the full PAM classification screen for a given protein/step.

        This replicates the computational screening that underlies the experimental
        data shown in Figure 3. It uses the trained model to score/evaluate all
        possible PAM sequences.

        Returns
        -------
        pandas.DataFrame: Complete classification table with predictions
                              (includes a 'Y_class' column marking '+' for cleaved).
        """
        model = self.load_model()
        df, _ = swiper.classify_all(model, n=self.n)
        return(df)

    def consensus(self, fnmask):
        """
        The main method for generating a consensus.

        This method executes the full workflow in one call:
        1. Screens all PAMs.
        2. Isolates only the "cleaved" (positive) predictions.
        3. Derives the degenerate consensus motif.
        4. Draws the sequence logo and saves it as a PNG.

        Parameters
        ----------
        fnmask (str): Base directory or filename prefix for the output PNG.
                          The final path will be: {fnmask}+.logo.png

        Returns
        -------
        consensus (str): The degenerate motif shown in the logo.
        """
        df = self.screen()
        cleaved = df[df["Y_class"] == "+"]
        consensus = convergence.get_consensus(cleaved["PAM"].values)
        ModelResult.draw_logo(cleaved, fnmask+".logo.png")
        return(consensus)
    
    
def get_consensus(TMP, SIZE, N_VALUE, what, file):
    """
    Prepares a temporary directory, extracts a zip file, identifies the final model file, and computes the consensus sequence.

    Parameters
    ----------
    TMP (str): The base path for the temporary directory. 
    SIZE (int): The model size. 
    N_VALUE (float): The N-value parameter.
    what (str): The specific model identifier. 
    file (str): The path to the file containing the model data. Performs a convergence analysis using the provided history data. Generates violation violin plots and saves them to the specified temporary directory. Calculates and returns the computed p-values.

:param TMP: The base path for the temporary directory where the plot is saved. :param history: The history data to be analyzed. :return: The computed p-values.
    
    Returns
    -------
    The consensus sequence.
    """
    os.system("mkdir "+TMP)
    os.system("mkdir "+TMP+"/"+what)
    os.system("unzip "+file+" -d "+TMP+"/"+what)
    models = [a for a in os.walk(TMP+"/"+what)][0][2]
    models = list(filter(lambda x: ".pkl" in x, models))
    final = get_final(models)
    print(final)
    path = TMP+"/"+what+"/"+final
    mr = ModelResult(path, n=SIZE, n_value=N_VALUE)
    consensus = mr.consensus(TMP+"/"+what)
    return(consensus)

def get_convergence(TMP, history):
    """
    Performs a convergence analysis using the provided history data. 
    Generates violation violin plots and saves them to the specified temporary directory. 
    Calculates and returns the computed p-values.
    
    Parameters
    ----------
    TMP (str): The base path for the temporary directory where the plot is saved. 
    history (np.array): The history data to be analyzed. 
    
    Returns
    -------
    The computed p-values.
    """
    lefts = convergence.generate_H0(10000, lags=[1,2,3,4,5], alphabet=convergence.IUPAC_DNA_no_gap)
    rights = convergence.make_rights(lefts, alphabet=convergence.IUPAC_DNA_no_gap)
    history_stat, history_H0, bvs, st, sigs = convergence.run_test(
        history, (lefts, rights), convergence.statistic, pval_threshold=0.05
    )
    plt.rcParams["figure.figsize"] = 10,10
    plt.rcParams["figure.dpi"] = 300
    fig, ax = plt.subplots()
    pvalues = np.array([convergence.compute_pvalue(a, history_H0) for a in history_stat])
    convergence.plot_test_violins(history, history_stat, history_H0, bvs, ax, 0.05, location="best", rotation=0)
    plt.savefig(TMP+"/convergence.png")
    return(pvalues)

def get_final(L):
    """
    Return the element from list L with the highest numeric suffix.

    Parameters
    ----------
    L (list): A list of strings.

    Returns
    -------
    str: The string from L with the highest integer suffix. 
    """
    return max(L, key=lambda s: int(op.splitext(s)[0].split('_')[1]))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Visualize swiper output: plot consensi and logos from left.zip and (optional) right.zip and the convergence statistics to output.pdf'
    )
    parser.add_argument(
        "--name",
        type=str,
        required=True,
        help="The experiment name"
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="The output filename"
    )
    parser.add_argument(
        "--experiments",
        type=str,
        required=True,
        help="exp1.zip;exp2.zip;exp3.zip"
    )
    parser.add_argument(
        "--size",
        type=int,
        default=7,
        help="The size parameter for the consensus generation"
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=2347,
        help="The seed for PRNG"
    )
    parser.add_argument(
        "--n_value",
        type=float,
        default=0.1,
        help="The n_value parameter for the consensus generation"
    )
    args = parser.parse_args()
    SIZE = int(args.size)
    N_VALUE = float(args.n_value)
    np.random.seed(args.seed)
    rnd = "".join(np.random.choice(np.arange(10), size=10).astype(str))
    TMP = "/tmp/swiper_"+rnd
    experiments = args.experiments.split(",")
    consensi = []
    for i,experiment in enumerate(experiments):
        current_consensus = get_consensus(TMP, SIZE, N_VALUE, str(i+1), experiment)
        current = (current_consensus, TMP+"/"+str(i+1)+".logo.png")
        consensi.append(current)
    convergence_ = None
    if len(consensi) > 1:
        pval = get_convergence(TMP, np.array([a[0] for a in consensi]).reshape(-1,1))
        convergence_ = (str(pval), TMP+"/convergence.png")
    o = Output(
        args.name, consensi, convergence_
    )
    o_md = o.out()
    output_fn = args.output
    with open(output_fn+".md", "w") as oh:
        oh.write(o_md)
    os.system("pandoc "+output_fn+".md -o "+output_fn+".pdf")