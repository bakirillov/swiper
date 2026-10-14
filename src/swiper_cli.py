"""Swiper CLI -- an Acive Learning application to guide PAM discovery experiments.
This file provides the code to run training and provide suggestions.
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
from weblogolib import *
from functools import partial
import matplotlib.pyplot as plt
from io import BytesIO, StringIO
from modAL.models import ActiveLearner
from sklearn.metrics import matthews_corrcoef
from modAL.batch import uncertainty_batch_sampling
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.metrics import accuracy_score, balanced_accuracy_score



def runs_by_what(fns, y_cols):
    """
    Group a collection of function (or file) names by the target column they
    refer to.

    Parameters
    ----------
    fns : Iterable[str]
        Iterable containing the names of the files to be
        inspected.
    y_cols : Iterable[str]
        Iterable of target column names (for example ``["AUC", "TPR"]``).

    Returns
    -------
    dict[str, list[str]]
        A mapping where the key is a column name from ``y_cols`` and the
        value is a list of every entry in ``fns`` that contains that
        column name somewhere in its string representation.
    """
    r = {a: [] for a in y_cols}
    for a in r:
        r[a] = list(filter(lambda x: str(a) in x, fns))
    return r


def compute_steps(v, n=50):
    """
    Generate a sequence of evenly‑spaced "step" values within the range
    spanned by the input array.

    Parameters
    ----------
    v : array‑like, shape (m,)
        Numerical data from which the minimum and maximum are
        derived.  The function does not sort the input; only the
        extreme values matter.
    n : int, optional (default: 50)
        Desired number of points *insid* the interval.  The actual
        number returned is `n‑1` because the first point is discarded.

    Returns
    -------
    numpy.ndarray, shape (n‑1,)
        A linearly spaced array that starts just after the minimum
        of `v` and ends just before the maximum of `v`.  The spacing
        between successive elements is `(max(v) - min(v)) / n`.
    """
    # Determine the numeric bounds of the input data.
    start = np.min(v)
    stop = np.max(v)
    # Create `n` points in the half‑open interval [start, stop).
    # The `endpoint=False` flag ensures that `stop` itself is not
    # included in the array.
    raw_steps = np.linspace(start, stop, n, endpoint=False)
    # Drop the first element to avoid
    # using the exact lower bound in subsequent calculations.
    return raw_steps[1:]


def random_query_strategy(model, X, n_instances):
    """
    Randomly select a batch of unlabeled instances for labeling.

    In an active‑learning pipeline the query strategy decides which
    examples from the unlabeled pool should be sent to an oracle (e.g.
    a lab experiment) next.  This implementation simply samples
    uniformly at random, which is useful as a baseline.

    Parameters
    ----------
    model : object
        The trained classifier (unused in this strategy but kept for
        API compatibility with more sophisticated strategies that
        might inspect the model).
    X : numpy.ndarray, shape (n_samples, n_features)
        Feature matrix of the unlabeled instances.
    n_instances : int
        Number of examples to query in this batch.

    Returns
    -------
    tuple
        * ``query_idx`` : 1‑D integer array of shape (n_instances,)
          containing the indices of the selected instances in ``X``.
        * ``X[query_idx]`` : the corresponding feature vectors.
    """
    # Randomly choose *n_instances* indices from the pool.
    # ``replace=True`` permits duplicate selections when the pool is
    # smaller than ``n_instances``.
    query_idx = np.random.choice(np.arange(X.shape[0]), replace=True, size=n_instances)
    # Return the indices and the corresponding feature vectors.
    return query_idx, X[query_idx]


def make_new_sheet(Ns, PAMs):
    """
    Create a fresh experiment sheet (pandas DataFrame) from raw lists.

    In the context of the project each experiment is stored as a sheet
    in an Excel workbook.  A sheet contains three columns:

        1. ``#``      – the sequence index (e.g. 1, 2, …)
        2. ``PAM``    – the nucleotide string that was tested
        3. ``Outcome`` – the experimental result (initially unknown)

    Parameters
    ----------
    Ns : array‑like
        Sequence numbers (typically 0‑based or 1‑based integers).
    PAMs : array‑like
        Corresponding PAM strings of equal length to ``Ns``.

    Returns
    -------
    pandas.DataFrame
        A new dataframe that can be written as a new sheet in the
        historical results workbook.  All entries in the ``Outcome``
        column are initialized to the empty string, signalling that
        the experiments have not yet been performed.
    """
    # ------------------------------------------------------------------
    # 1. Sanity check – a one‑to‑one mapping between indices and PAMs
    # ------------------------------------------------------------------
    assert len(Ns) == len(PAMs)
    # ------------------------------------------------------------------
    # 2. Build the dataframe with the required three columns.
    #    * ``Outcome`` is filled with an empty string for every row.
    # ------------------------------------------------------------------
    return (
        pd.DataFrame(
            {
                "#": Ns,          # sequence identifier
                "PAM": PAMs,      # tested sequence
                "Outcome": [""] * len(Ns)  # placeholder for future results
            }
        )
    )


def make_history(history, temp_folder):
    """
    Persist a collection of data‑frames to an Excel workbook.

    The helper is used by the experiment framework to snapshot the
    history of all runs (e.g. performance metrics, parameter
    combinations, etc.) in a single Excel file that can later be
    inspected manually or read back into Python.

    Parameters
    ----------
    history : dict[str, pandas.DataFrame]
        Mapping from a logical sheet name to a DataFrame that will be
        written to that sheet.  The key string is used verbatim as the
        Excel worksheet name.
    temp_folder : str or pathlib.Path
        Directory where the temporary workbook should be created.
        The function does not create the folder; it assumes that
        `temp_folder` already exists.

    Returns
    -------
    None
        The function writes the file in place and closes the writer
        automatically.  No object is returned.
    """
    # Build an ``ExcelWriter`` that will write to
    #   <temp_folder>/result.xlsx
    writer = pd.ExcelWriter(op.join(temp_folder, "result.xlsx"),
                            engine="xlsxwriter")
    # Iterate over the keys in the `history` dictionary.
    # Each key is treated as the name of an Excel worksheet.
    # The associated DataFrame is written into that sheet.
    for a in history:
        history[a].to_excel(writer, sheet_name=a)
    # Close the writer so that the file is flushed and the workbook
    # is finalized.  No return value is needed.
    writer.close()


def classify_X(classifier, X, df):
    """
    Apply a trained model to a new set of features and return two dataframes:

    1. A dataframe that contains all the original rows plus a column
       ``Y_class`` with the class prediction for each sample.

    2. A subset dataframe containing only the samples that the
       classifier has flagged as most informative for the next
       active‑learning query (the “query set”).

    Parameters
    ----------
    classifier : sklearn‑like estimator with a ``predict`` method and
                 an optional ``query`` method
        The model that was previously trained on a training set.  It is
        expected to implement the active learning API used in this
        project, i.e.:

        * ``predict(X)`` returns a 1‑D array of class labels.
        * ``query(X)`` returns a tuple ``(idx, instance)`` where
          * ``idx``  – indices of the most uncertain samples in ``X``
          * ``instance`` – the corresponding feature vectors (not used
          here but returned for consistency).

    X : np.ndarray, shape (n_samples, n_features)
        Feature matrix for the new data to classify.

    df : pandas.DataFrame
        A dataframe that must contain the same number of rows as
        ``X`` (in the same order).  It is usually created by reading
        a FASTA or similar file.  Only the column names of ``df`` are
        preserved – new columns are added in‑place.

    Returns
    -------
    tuple
        * ``df`` (modified): the original dataframe with an added
          column ``Y_class`` containing the predicted class for each
          row.
        * ``df.iloc[query_idx]``: a dataframe with the rows that were
          selected by ``classifier.query``.  These rows are the
          most uncertain samples and will typically be sent back to
          the experimenter for labeling.

    Example
    -------
    >>> df_pool = pd.DataFrame({'PAM': ['ACGT', 'CGTA', 'GTAC']})
    >>> X_pool = np.stack([onehot(s) for s in df_pool['PAM']])
    >>> model = make_RF()
    >>> model = train(X_train, y_train, model, uncertainty_batch_sampling, 10, tmp)
    >>> df_full, df_query = classify_X(model, X_pool, df_pool)
    >>> df_full.head()
    """
    # Make a local copy of the dataframe so that we do not modify the
    # original object that was passed in (the slicing syntax ``df[:]``
    # creates a shallow copy; pandas may return a view or a copy depending
    # on the underlying data layout).
    df = df[:]
    # 1. Predict the class for every example in X
    y_full_class = classifier.predict(X)
    # 2. Attach those predictions as a new column
    df["Y_class"] = y_full_class
    # 3. Ask the classifier which samples should be queried next.
    #    The ``query`` method typically returns indices of the
    #    most uncertain samples and their corresponding feature
    #    vectors.  Only the indices are needed here.
    query_idx, query_instance = classifier.query(X)
    # 4. Return the enriched dataframe and a subset containing the
    #    query samples.  ``df.iloc[query_idx]`` selects the rows by
    #    integer position.
    return (df, df.iloc[query_idx])


def classify_all(classifier, n=7):
    """
    Classify every possible PAM sequence of a given length.

    The function is useful when the experimenter wants to know how the
    current model would classify all sequences in the combinatorial
    space before deciding which ones to test next.

    Parameters
    ----------
    classifier : sklearn‑like estimator
        A model that implements ``predict`` and ``query`` (the
        active‑learning API).  It is expected to have been trained
        previously on a labelled dataset.

    n : int, optional (default: 7)
        Length of the PAM sequences to enumerate.  The total number
        of combinations is ``4**n`` (since each position can be
        A/T/G/C).

    Returns
    -------
    tuple
        * A pandas.DataFrame containing a column ``PAM`` with all
          sequences and an added column ``Y_class`` with the model’s
          predictions (returned by :func:`classify_X`).
        * A second DataFrame that contains only the rows that the
          classifier has selected as the most informative queries
          (the query set).
    """
    # ------------------------------------------------------------------
    # 1. Enumerate all possible PAM sequences of length ``n``.
    # ------------------------------------------------------------------
   # print("Sequence size " + str(n))
    # Create a list that contains the same 4‑letter alphabet repeated
    # ``n`` times – e.g. ["ATGC", "ATGC", ..., "ATGC"] (n times).
    # ``itertools.product`` will then produce the Cartesian product
    # over this list, giving every combination of A,T,G,C at each
    # position.
    base = ["ATGC"] * n
    # ``product(*base)`` yields tuples like ('A', 'T', 'G'), which are
    # concatenated into a single string: "ATG".
    pam_strings = ["".join(a) for a in itertools.product(*base)]
    # ------------------------------------------------------------------
    # 2. Put the sequences into a DataFrame so that the original
    #    column names are preserved and we can add predictions later.
    df = pd.DataFrame({"PAM": np.array(pam_strings)})
    # ------------------------------------------------------------------
    # 3. Convert the strings into the numeric feature matrix expected
    #    by the model.  ``onehot`` is assumed to produce a fixed‑size
    #    vector for each sequence.
    X = np.stack(df["PAM"].apply(onehot).values).astype(np.float32)
    # ------------------------------------------------------------------
    # 4. Delegate the actual prediction + active‑learning query to
    #    ``classify_X`` so that the same logic is applied.
    return classify_X(classifier, X, df)


def save_fasta(fasta_fn, sequences):
    """
    Write a list of sequences to a FASTA file.

    Parameters
    ----------
    fasta_fn : str
        Path (including filename) where the FASTA file should be written.
    sequences : Iterable[str]
        A list or other iterable of raw DNA/RNA/protein sequences.
        Each element must be a string containing only the sequence
        characters (e.g. "ACGT", "ATGCA").

    Returns
    -------
    None
        The function writes the file and closes the handle
        without returning anything.
    """
    # Open the target file for writing (overwrites if it already exists).
    oh = open(fasta_fn, "w")
    # Enumerate over the supplied sequences.  ``i`` becomes the
    #   header identifier, ``a`` is the actual sequence.
    for i, a in enumerate(sequences):
        oh.write(">" + str(i) + "\n" + a + "\n")
    # Explicitly close the file handle.
    oh.close()


def make_logo(fasta_fn, png_fn, files=True):
    """
    Generate a sequence logo image from a FASTA file and write it to disk.
    
    Parameters
    ----------
    fasta_fn : str | file‑like
        Path to a FASTA file or an already opened file object.
        If ``files`` is ``False`` the argument is expected to be a
        file‑like object (i.e. something that implements ``read``).
    png_fn : str
        Path where the resulting PNG image will be written when
        ``files`` is ``True``.  Ignored otherwise.
    files : bool, optional
        * ``True`` – treat ``fasta_fn`` as a filename and write the
          PNG to ``png_fn``.  Return ``None``.
        * ``False`` – treat ``fasta_fn`` as an open file‑like object
          and return the PNG byte stream directly.
    
    Returns
    -------
    None | bytes
        ``None`` if the image was written to disk, otherwise a
        ``bytes`` object containing the PNG data.
    """
    # ------------------------------------------------------------------
    # Step 1 – open the FASTA input (either a file path or a file‑like).
    # ------------------------------------------------------------------
    if files:
        # ``fasta_fn`` is a path string: open it for reading.
        fin = open(fasta_fn)
    else:
        # Assume the caller supplied an already opened file‑like object.
        fin = fasta_fn
    # ------------------------------------------------------------------
    # Step 2 – parse the FASTA content into a list of raw sequences.
    # ------------------------------------------------------------------
    seqs = read_seq_data(fin)
    # ------------------------------------------------------------------
    # Step 3 – build the data structures that describe the logo.
    # ------------------------------------------------------------------
    logodata = LogoData.from_seqs(seqs)          # per‑position frequencies
    logooptions = LogoOptions()                  # default visual options
    logoformat = LogoFormat(logodata, logooptions)  # how to render
    # ------------------------------------------------------------------
    # Step 4 – render the logo as a PNG byte stream.
    # ------------------------------------------------------------------
    png = png_formatter(logodata, logoformat)
    # ------------------------------------------------------------------
    # Step 5 – write to disk or return the PNG bytes.
    # ------------------------------------------------------------------
    if files:
        # Write the binary PNG data to ``png_fn``.
        with open(png_fn, "wb") as oh:
            oh.write(png)
        return None
    else:
        # Caller wants the PNG data as a bytes object.
        return png


def toheno(u):
    """
    Decode a flattened one‑hot array back into a nucleotide string.

    The project encodes the sequences as a 1‑D array of length
    ``4 * L`` where each group of four consecutive elements represents
    a single base in a one‑hot scheme.  ``toheno`` reverses that process
    and yields the corresponding string of bases.

    Parameters
    ----------
    u : numpy.ndarray
        A 1‑D NumPy array containing a concatenated one‑hot encoding.
        The array must have a length that is a multiple of four.
        Example shapes:

        * ``(4*L,)`` – a single sequence flattened into one dimension
        * ``(n,)``   – same as above; the function does not care how
          the array was produced, as long as it satisfies the shape
          requirement.

        Each 4‑element block must match one of the four canonical
        encodings.  If an element does not correspond to a known
        base, a ``KeyError`` will be raised.

    Returns
    -------
    str
        The decoded sequence, e.g. ``array([1,0,0,0, 0,0,0,1])`` →
        ``"AT"``.
    """
    # ------------------------------------------------------------------
    # Mapping from the string representation of a 4‑element one‑hot
    # vector to the corresponding nucleotide.
    # ------------------------------------------------------------------
    encoding = {
        "1,0,0,0": "A",
        "0,0,0,1": "T",
        "0,1,0,0": "G",
        "0,0,1,0": "C"
    }
    # ------------------------------------------------------------------
    # Reshape `u` from a flat vector into an (L, 4) matrix.
    #   u.shape[0]   → 4*L
    #   int(u.shape[0] / 4) → L
    # The `order='f'` (Fortran style) ensures that the first 4
    # elements of the flat array form the first row of the matrix,
    # matching the original packing of the one‑hot representation.
    # We then transpose (`.T`) so that each row corresponds to one
    # base.  After the transpose, `h` has shape (L, 4).
    # ------------------------------------------------------------------
    h = u.reshape(4, int(u.shape[0] / 4), order="f").T
    # ------------------------------------------------------------------
    # Convert each row of the matrix into a comma‑separated string
    # (e.g. [1,0,0,0] → "1,0,0,0") so it can be used as a dictionary key.
    # ``int(a)`` guarantees we use 0/1 values (in case `u` is float).
    # ------------------------------------------------------------------
    e = [",".join([str(int(a)) for a in c]) for c in h]
    # ------------------------------------------------------------------
    # Translate each encoded string back into its nucleotide and
    # concatenate them to form the final sequence.
    # ------------------------------------------------------------------
    r = "".join([encoding[a] for a in e])
    return r


def onehot(u):
    """
    Convert a sequence of nucleotides (or integer codes) into a
    concatenated one‑hot encoding.

    Parameters
    ----------
    u : iterable
        An iterable of characters or integer codes representing a
        DNA/RNA sequence.  Accepted values are:

        * Integers 0–4:
            0 → [0, 0, 0, 0]   (unknown / gap)
            1 → [1, 0, 0, 0]   (A)
            2 → [0, 0, 0, 1]   (T)
            3 → [0, 1, 0, 0]   (G)
            4 → [0, 0, 1, 0]   (C)
        * Single‑character strings, case‑insensitive:
            'A','a' → [1,0,0,0]
            'T','t' → [0,0,0,1]
            'G','g' → [0,1,0,0]
            'C','c' → [0,0,1,0]
            'N','n' → [0,0,0,0]   (ambiguous / unknown)

        Any value not in the table will raise a ``KeyError``.

    Returns
    -------
    numpy.ndarray
        A one‑dimensional NumPy array containing the flattened
        one‑hot vectors for each element of u, concatenated in
        the original order.  For a sequence of length L, the
        returned shape is ``(L*4,)``; e.g. ``['A','T']`` → ``[1,0,0,0,0,0,0,1]``.
    """
    # ------------------------------------------------------------------
    # Define the mapping from bases / integer codes to 4‑dimensional
    # one‑hot vectors.  The dictionary contains all allowed keys.
    # ------------------------------------------------------------------
    encoding = {
        1: [1,0,0,0],
        2: [0,0,0,1],
        3: [0,1,0,0],
        4: [0,0,1,0],
        0: [0,0,0,0],
        "A": [1,0,0,0],
        "T": [0,0,0,1],
        "G": [0,1,0,0],
        "C": [0,0,1,0],
        "N": [0,0,0,0],
        "a": [1,0,0,0],
        "t": [0,0,0,1],
        "g": [0,1,0,0],
        "c": [0,0,1,0],
        "n": [0,0,0,0]
    }
    # ------------------------------------------------------------------
    # Translate each element of u using the mapping.  The list
    # comprehension yields a list of 4‑element vectors.
    # ``sum(..., [])`` concatenates all those lists into a single
    # flat list; it is a quick trick equivalent to ``itertools.chain``.
    # ------------------------------------------------------------------
    r = np.array(sum([encoding[a] for a in u], []))
    # ------------------------------------------------------------------
    # Return the flattened NumPy array.
    # ------------------------------------------------------------------
    return r


def make_RF():
    """
    Create and return a fresh scikit‑learn ``RandomForestClassifier``.

    Returns
    -------
    sklearn.ensemble.RandomForestClassifier
        An un‑trained random‑forest model ready for ``fit``/``predict``.

    Notes
    -----
    * No arguments are passed to the constructor, so all parameters
      assume their default values.  For reference the defaults are
      roughly:
         RandomForestClassifier(
             n_estimators=100,          # number of decision trees
             criterion='gini',          # impurity measure
             max_depth=None,            # grow trees fully
             min_samples_split=2,       # node split threshold
             min_samples_leaf=1,        # leaf node size
             max_features='auto',       # features considered at each split
             bootstrap=True,            # bootstrap samples
             n_jobs=None,               # parallelism (uses all CPUs by default)
             random_state=None          # reproducibility seed
         )
    """
    # Instantiate a RandomForestClassifier with default hyper‑parameters.
    rf = RandomForestClassifier()
    return rf


def train(X, Y, model_fun, query_fun, query_size, temp_folder):
    """
    Train an Active Learner on a collection of data splits.

    Parameters
    ----------
    X : dict[str, array-like]
        A mapping from a chunk identifier to the input features.
        Each entry must be a 2‑D array (samples, features).  The
        special key ``"1"`` denotes the initial training set.
    Y : dict[str, array-like]
        Corresponding labels for the data in `X`.  The shape of
        `Y[key]` must match the number of rows in `X[key]`.
    model_fun : callable
        A no‑argument function that returns a fresh estimator
        (e.g. `lambda: RandomForestClassifier()`).
    query_fun : callable
        A query strategy function that takes a model and the
        training data and returns a list of indices of samples
        to label.  It must accept the keyword argument
        ``n_instances`` to control how many instances to query.
    query_size : int
        The number of instances to request from the query strategy
        during the initialisation phase.
    temp_folder : str
        Path to a directory where intermediate pickles will be
        written.  The folder must exist and be writable.

    Returns
    -------
    ActiveLearner
        The fully trained learner after all chunks have been
        processed.  The returned object is the same one that is
        pickled after each iteration.
    """
    # --------------------------------------------------------------------
    # Print the size of the query used for the first round
    # --------------------------------------------------------------------
    print("Query size " + str(query_size))
    # --------------------------------------------------------------------
    # Separate the *initial* training set (key "1") from the rest
    # --------------------------------------------------------------------
    X_training = X["1"]          # features for the first chunk
    y_training = Y["1"]          # corresponding labels
    # --------------------------------------------------------------------
    # Build a preset query strategy that always asks for
    #      `query_size` instances.
    # --------------------------------------------------------------------
    # `partial` freezes the `n_instances` argument so that the
    # ActiveLearner can call it with just the model and data.
    preset_batch = partial(query_fun, n_instances=query_size)
    # --------------------------------------------------------------------
    # Initialise the ActiveLearner
    # --------------------------------------------------------------------
    #  * `estimator`   – fresh model created by `model_fun()`
    #  * `X_training`/`y_training` – the initial labelled set
    #  * `query_strategy` – the preset batch query strategy
    classifier = ActiveLearner(
        estimator=model_fun(),
        X_training=X_training,
        y_training=y_training,
        query_strategy=preset_batch
    )
    # --------------------------------------------------------------------
    # Persist the initial learner state
    # --------------------------------------------------------------------
    with open(op.join(temp_folder, "model_1.pkl"), "wb") as oh:
        pkl.dump(classifier, oh)
    print("Processing 1")
    # --------------------------------------------------------------------
    # Iteratively teach the learner with each remaining chunk
    # --------------------------------------------------------------------
    for a in X.keys():
        if a != "1":                     # skip the initial chunk
            print("Processing " + a)
            classifier.teach(X[a], Y[a])  # add new labelled data
            # Save the updated model after this teach step
            with open(op.join(temp_folder, f"model_{a}.pkl"), "wb") as oh:
                pkl.dump(classifier, oh)
    # --------------------------------------------------------------------
    # Return the fully trained classifier
    # --------------------------------------------------------------------
    return classifier


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d", "--dataset",
        dest="dataset",
        action="store", 
        help="set the dataset with previous experiment results", 
        default="initial.xlsx",
    )
    parser.add_argument(
        "-p", "--pool",
        dest="pool",
        action="store", 
        help="set the pool to select data from", 
        default="initial.fasta",
    )
    parser.add_argument(
        "-s", "--seed",
        dest="seed",
        action="store", 
        help="set the seed for random number generation", 
        default="10",
    )
    parser.add_argument(
        "-q", "--query-size",
        dest="query_size",
        action="store", 
        help="set the query size",
        default="18"
    )
    parser.add_argument(
        "-o", "--output",
        dest="output",
        action="store", 
        help="set the output data filename", 
        default="output.zip",
    )
    parser.add_argument(
        "-r", "--random",
        dest="random",
        action="store_true", 
        help="use random sampling",
    )
    parser.add_argument(
        "-a", "--archive",
        dest="archive",
        action="store_true", 
        help="make an archive",
    )
    parser.add_argument(
        "-t", "--temp-dir",
        dest="tempdir",
        action="store", 
        help="temporary directory",
    )
    args = parser.parse_args()
    seed = int(args.seed)
    np.random.choice(seed)
    xls = pd.ExcelFile(args.dataset)
    history = pd.read_excel(xls, None)
    history = {a: history[a][["#", "PAM", "Outcome"]] for a in history}
    query_size = int(args.query_size)
    X = {a: np.stack([onehot(b).astype(np.float32) for b in history[a]["PAM"]]) for a in history}
    Y = {a: history[a]["Outcome"].values for a in history}
    temp_folder = args.tempdir
    os.makedirs(temp_folder, exist_ok=True)
    query_fun = uncertainty_batch_sampling if not args.random else random_query_strategy
    print(args.random)
    if args.random:
        print("Using random sampling")
    else:
        print("Using active learning")
    model_fun = make_RF
    print("X,Y", [X[a].shape for a in X], [Y[a].shape for a in Y])
    model = train(X, Y, model_fun, query_fun, query_size, temp_folder)
    next_sheet_num = str(int(np.max([int(a) for a in list(history.keys())])+1))
    if args.pool == "None":
        print("Using entire pool")
        full_data, all_df_q = classify_all(model, len(history["1"]["PAM"][0]))
        new_sheet = make_new_sheet(all_df_q.index.values, all_df_q["PAM"].values)
        history[next_sheet_num] = new_sheet
    else:
        print("Using "+args.pool)
        pool = []
        with open(args.pool, "r") as ih:
            for a in ih.readlines():
                if a[0] != ">":
                    pool.append(a.replace("\n", ""))
        pool = np.array(pool)
        pooldf = pd.DataFrame({"PAM": pool})
        poolX = np.stack([onehot(a) for a in pool]).astype(np.float32)
        full_data, all_df_q = classify_X(model, poolX, pooldf)
        new_sheet = make_new_sheet(all_df_q.index.values, all_df_q["PAM"].values)
        history[next_sheet_num] = new_sheet
    make_history(history, temp_folder)
    if args.archive:
        shutil.make_archive(args.output, 'zip', temp_folder)