"""This file provides the code to compute statistical tests for active learning convergence.
This file is to be imported as a library and provides the following functions:

    get_consensus: Converts a list of PAM sequences into a degenerate IUPAC consensus string 
        to visualize the cleavage preference.
    consensus_score: 
        Computes a similarity score between a reference (real) PAM sequence and a
        predicted degenerate consensus sequence.
    lead_trail:
        Splits the PAM candidate into leading and trailing part.
    get_probs:
        Computes probability distribution of a base selection over the alphabet 
        letters given the counts and the leading part of the PAM candidate.
    generate_H0
        Generate left samples from H0 (the active learning did not converge) for the statistical test.
    mutate:
        Return s with randomly changed letters.
    make_rights:
        Generate right strings for H0 sample.
    statistic:
        Statistics according to the Eq. 2 from the paper.
    get_boundary:
        Compute boundary values for the given H0 sample.
    compute_pvalue:
        Compute p-value for a given statistics value.
    run_test:
        Runs statistical test for given experimental history.
    plot_test_violins:
        Violin plots for difference between experimental and H0 samples.
"""



import itertools
import numpy as np
from Bio import motifs
from Bio.Seq import Seq
from collections import Counter
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


IUPAC_DNA_no_gap = "ACGTRYSWKMBDHVN"


def get_consensus(l):
    """
    Converts a list of PAM sequences into a degenerate IUPAC consensus string 
    to visualize the cleavage preference.

    Parameters
    ----------
        l (list of str): List of DNA sequences (the "cleaved" PAMs from the screen).

    Returns
    -------
        str: Degenerate consensus sequence (e.g., "NGG").
    """
    consensus = str(motifs.create([Seq(a) for a in l]).degenerate_consensus)
    return(consensus)


def consensus_score(real, predicted, n_value=0.1):
    """
    Computes a similarity score between a reference (real) PAM sequence and a
    predicted degenerate consensus sequence.

    The scoring is position‑wise and accounts for degenerate IUPAC bases:
    - Positions where the reference is 'N' (any base) are considered less informative
      and are scored with a reduced weight (n_value) unless they occur at the
      beginning of the sequence (left‑trailing Ns), which receive full score (1).
    - Exact matches to non‑'N' bases yield a score of 1.
    - Mismatches yield 0.

    This scoring scheme biases toward exact conservation of known specific bases,
    while giving some credit for matches at positions that are known to be flexible
    (N) in the reference. The trailing Ns after the first specific base are
    down‑weighted because they represent less constrained positions.

    Parameters
    ----------
        real (str): The reference/wild‑type PAM sequence (e.g., "NGG").
        predicted (str): The degenerate consensus sequence derived from cleavage data
                         (e.g., "RGG").
        n_value (float, optional): Weight assigned to matches where the reference
                                   base is 'N' and the position is not left‑trailing.
                                   Defaults to 0.1.

    Returns
    -------
        float: The cumulative similarity score. Higher values indicate closer
               agreement with the reference, with a maximum possible score equal to
               the length of the sequence (if all positions match perfectly with
               weight 1).
    """
    score = 0
    left_trail = True
    for a,b in zip(real, predicted):
        add = 0
        if a != "N":
            left_trail = False
        if a == b:
            if a == "N":
                if left_trail:
                    add = 1
                else:
                    add = n_value
            else:
                add = 1
        score += add
    return(score)


def lead_trail(template):
    """
    Splits the PAM candidate into leading and trailing part.

    Parameters
    ----------
        template (string): PAM candidate to split

    Returns
    -------
        (leading_part, trailing_part) (Tuple of strings): leading part and the Ns at the end of the PAM
    """
    # Find the lagging part (trailing Ns only)
    # We need to find the last non-N character, then everything after it is trailing Ns
    last_non_n_index = -1
    for i in range(len(template) - 1, -1, -1):
        if template[i] != 'N':
            last_non_n_index = i
            break
    # Split template into leading and trailing parts
    leading_part = template[:last_non_n_index + 1]  # Up to and including the last non-N
    trailing_part = template[last_non_n_index + 1:]  # All trailing Ns after the last non-N
    return(leading_part, trailing_part)


def get_probs(alphabet_letters, letter_counts, leading_part):
    """
    Computes probability distribution of a base selection over the alphabet 
    letters given the counts and the leading part of the PAM candidate.

    Parameters
    ----------
        alphabet_letters(Iterable): the letters in the alphabet
        letter_counts([int]): counts of each letter in the template
        leading part(string): the leading part of the PAM candidate
    Returns
    -------
        [float]: probabilities to be picked for each letter
    """
    probs = []
    # Small probability for letters not in template
    uniform_prob = 0.01 / len(alphabet_letters)
    # Calculate total probability mass for template letters
    num_unique_in_template = len(letter_counts)
    num_not_in_template = len(alphabet_letters) - num_unique_in_template
    # Adjust template_prob_mass to account for uniform_prob
    template_prob_mass = max(0, 1.0 - (uniform_prob * num_not_in_template))
    for letter in alphabet_letters:
        if letter in letter_counts:
            # Scale template frequency to fit remaining probability mass
            prob = template_prob_mass * (letter_counts[letter] / len(leading_part))
        else:
            prob = uniform_prob
        probs.append(prob)
    # Normalize to ensure sum = 1
    total = sum(probs)
    if total > 0:
        probs = [p / total for p in probs]
    else:
        # Fallback to uniform distribution if something went wrong
        probs = [1/len(alphabet_letters)] * len(alphabet_letters)
    return(probs)


def generate_H0(sample_size, alphabet, sequence_length=7, lags=[1,2,3,4], seed=None):
    """
    Generate left samples from H0 (the active learning did not converge) for the statistical test.
    
    We generate “left-hand” strings by drawing random sequences from a fixed alphabet
    (IUPAC nucleotides with no gaps) while deliberately inserting a random number of
    ambiguous nucleotides (‘N’s) at the terminus (to simulate the trailing part of PAM).
    By sampling the trailing part size uniformly from the provided set of lengths, we
    emulate the diversity of possible PAM sequences.

    Parameters
    ----------
        sample_size(int): number of generated left strings
        alphabet(str): the degenerate IUPAC alphabet for nucleotides
        sequence_length(int): the length of each PAM candidate
        lags([int]): the sizes of random trailing part
        seed(int): the random number generator seed

    Returns
    -------
        numpy.ndarray of strings: the generated sample of length sample_size
    """
    # Set the random seed if provided
    if seed is not None:
        np.random.seed(seed)
    # Get unique alphabet letters
    alphabet_letters = sorted(set(alphabet))
    # Count occurrences of each letter in leading part (including any Ns in the middle)
    # letter_counts = Counter(leading_part)
    # Create probability distribution:
    probs = [1/len(alphabet_letters)]*len(alphabet_letters)
    # Create mapping from index to letter
    letter_dict = {i: alphabet_letters[i] for i in range(len(alphabet_letters))}
    # Generate random strings
    generated_strings = []
    for _ in range(sample_size):
        # Generate random indices according to the probability distribution
        random_trail = np.random.choice(lags)
        random_indices = np.random.choice(
            len(alphabet_letters), 
            size=sequence_length-random_trail, 
            p=probs
        )
        # Convert indices to letters
        leading_random = ''.join([letter_dict[idx] for idx in random_indices])
        # Add the trailing Ns part
        trailing_part = "N"*random_trail
        full_string = leading_random + trailing_part
        generated_strings.append(full_string)
    return np.array(generated_strings)


def mutate(s, alphabet):
    """
    Return s with randomly changed letters.
    
    The selected positions are sampled without replacement to avoid repeated changes at the same site, 
    ensuring that each mutation event modifies distinct locus. For each chosen position, the function
    constructs a list of alternative characters that excludes the current base, making sure
    that a genuine change occurs. 
    
    Parameters
    ----------
        s (str): Input string
        alphabet (str): String containing allowed letters
    
    Returns
    -------
        str: Mutated string with random changes
    """
    # Convert to list for mutable operations
    s_list = list(s)
    n = len(s_list)
    # Random number of positions to change (1 to n)
    num_changes = np.random.randint(1, n + 1)
    # Random positions to change (without replacement)
    positions = np.random.choice(n, size=num_changes, replace=False)
    for pos in positions:
        # Current character at position
        current_char = s_list[pos]
        # Available characters (excluding current one)
        available_chars = [c for c in alphabet if c != current_char]
        if available_chars:  # If there are alternatives
            # Randomly choose a different character
            new_char = np.random.choice(available_chars)
            s_list[pos] = new_char
    return(''.join(s_list))


def make_rights(lefts, alphabet):
    """
    Generate right strings for H0 sample.
    
    The “right-hand” strings are generated by replacing a randomly chosen subset 
    of characters with alternative symbols from the same alphabet.
    
    Parameters
    ----------
        lefts (Iterable of str): Input left strings
        alphabet (str): String containing allowed letters
    
    Returns
    -------
        numpy.ndarray: Array of mutated right strings with random changes
    """
    return(np.array([mutate(s, alphabet) for s in lefts]))


def statistic(left, right):
    """
    Statistics according to the Eq. 2 from the paper.
    
    1. We compute consensus with each other, but, since vanilla consensus scores often lead to 
        premature termination at suboptimal results:
    2. We divide consensus scores by the leading part length to navigate past local maxima and:
    3. We penalize convergence if the leading part is too short (shorter than two nucleotides).
    
    Parameters
    ----------
        left (str): left string from the pair (n-th experiment).
        right (str): right string from the pair (n+1-th experiment).
    
    Returns
    -------
        float: value of the statistic.
        
    """
    left_lead, _ = lead_trail(left)
    right_lead, _ = lead_trail(right)
    consensus_with_each_other = consensus_score(left, right)#+10**-5)
    max_len = np.max([len(right_lead), len(left_lead)])
    if max_len <= 2:
        max_len = 10**5
    return(consensus_with_each_other/max_len)

def get_boundary(history_H0, n_comparisons, pval_threshold=0.05):
    """Compute boundary values for the given H0 sample."""
    return(np.quantile(history_H0, 1-pval_threshold/n_comparisons))

def compute_pvalue(stat, H0):
    """Compute p-value for a given statistics value."""
    return(H0[H0 >= stat].shape[0]/H0.shape[0])

def run_test(history, LR, statistic_f, pval_threshold=0.05):
    """
    Runs statistical test for given experimental history.
    
    The test quantifies the similarity between two motifs adjusted by the quality of motif extraction, 
    and the p-value is computed with Monte Carlo simulation by estimation of probability for a random pair of PAM candidates 
    having the same or more extreme concordance under the null hypothesis that no true convergence has occurred.

    Parameters
    ----------
        history([pd.DataFrame]): experimental history for an ortholog.
        LR(([str], [str]): lefts and rights sampled from H0.
        statistics_f((str, str) -> float): function to compute statistics.
        pval_threshold(float): significance threshold.
    Returns
    ------
        history_stat(np.array(float)): statistics computed over experimental sample.
        history_H0(np.array(float)): statistics computed over H0 sample.
        boundary_values(np.array(float)): boundary values for the tests.
        sig_test(np.array(bool)): results of significance tests.
        significants(np.array((str, str))): only statistically significant pairs
    """
    lefts, rights = LR
    # Compute H0 statistics
    history_H0 = np.array([statistic_f(a,b) for a,b in zip(lefts, rights)])
    # Compute left-right pairs from the experiment and experimental statistics
    history_stat = []
    pairs = []
    for i in np.arange(history.shape[0]-1):
        history_stat.append(statistic_f(history.T[0][i], history.T[0][i+1]))
        pairs.append((i, i+1))
    pairs = np.array(pairs)
    history_stat = np.array(history_stat)
    # Compute boundary values
    boundary_values = np.array(
        [get_boundary(history_H0, n+1, pval_threshold) for n in np.arange(0,history_stat.shape[0])]
    )
    # Compute significance test
    sig_test = history_stat > boundary_values
    significants = pairs[np.where(sig_test)[0]]
    return(history_stat, history_H0, boundary_values, sig_test, significants)


def plot_test_violins(history, history_stat, history_H0, boundaries, ax, pval, location, rotation):
    """
    Violin plots for difference between experimental and H0 samples.

    Parameters
    ----------
        history([pd.DataFrame]): experimental history for an ortholog.
        history_stat(np.array(float)): statistics computed over experimental sample.
        history_H0(np.array(float)): statistics computed over H0 sample.
        boundary_values(np.array(float)): boundary values for the tests.
        ax(plt.Axis): Matplotlib axis.
        pval(float): p-value threshold.
        location(str): legend location on the plot.
        rotation(float): rotation of tick labels.
        

    Returns
    -------
        None
        
    """
    positions = np.arange(1, history_stat.shape[0]+1)
    observed_stats = []
    labels = []
    for i, c in enumerate(history_stat):
        observed_stats.append(c)
        current_pval = compute_pvalue(c, history_H0)
        labels.append(f"{i}: {history.T[0][i]}→{history.T[0][i+1]}\np={current_pval}")
    violin_parts = ax.violinplot([history_H0]*len(positions), positions=positions, vert=True, 
                             showmeans=False, showmedians=False, showextrema=False)
    for i, pc in enumerate(violin_parts['bodies']):
        pc.set_facecolor("black")
        pc.set_edgecolor('black')
        pc.set_edgecolors("black")
        pc.set_alpha(0.8)
    # Add observed statistic points
    for i, (pos, stat, bnd) in enumerate(zip(positions, observed_stats, boundaries)):
        ax.scatter(pos, stat, color="yellow", s=50, zorder=3, 
               marker="x", linewidth=2, label=f"Observed stat" if i == 0 else None)
        ax.scatter(pos, bnd, color="snow", s=50, zorder=3, 
               marker=".", linewidth=2, label=f"Boundary value at p="+str(pval) if i == 0 else None)
    ax.set_xlabel("Gel Pair Index")
    ax.set_ylabel("Statistic Value")
    #ax.set_title("Statistical test for convergence of active learning-derived PAM consensus")
    ax.set_xticks(positions)
    ax.set_xticklabels(labels, rotation=rotation, ha='center')
    ax.grid(True, alpha=0.3, linestyle='--', axis='y')
    legend_elements = [
        Line2D([0], [0], color='black', alpha=0.8, lw=4, label='H0: No Convergence'),
        Line2D([0], [0], marker='x', color='w', markerfacecolor='yellow', markersize=8, 
           markeredgecolor='yellow', markeredgewidth=2, label='Observed Statistic'),
        Line2D(
            [0], [0], marker=".", color="w", markerfacecolor='snow', alpha=0.8, lw=4, 
            label='Boundary value at p='+str(pval)
        ),
    ]
    ax.legend(handles=legend_elements, loc=location)