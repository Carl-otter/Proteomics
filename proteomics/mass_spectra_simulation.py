def calculate_mol_mass(peptide_seq, amino_acid_mass_dict):
    """
    takes a single peptide sequence and returns a dictionary `{peptide: mass}` by summing the amino-acid masses.

    Parameters
    ----------
    peptide_seq : str
        A string representing the peptide sequence.
    amino_acid_mass_dict : dict
        A dictionary mapping each amino acid to its corresponding mass in Daltons.
    
    Returns
    -------
    dict
        A dictionary containing the peptide sequence as the key and its calculated molecular mass as the value.
    """
    amino_acid_mass_dalton = {
    'A': 71.08, 'R': 156.19, 'N': 114.10, 'D': 115.09,
    'C': 103.15, 'E': 129.12, 'Q': 128.13, 'G': 57.05,
    'H': 137.14, 'I': 113.16, 'L': 113.16, 'K': 128.17,
    'M': 131.19, 'F': 147.18, 'P': 97.12, 'S': 87.08,
    'T': 101.11, 'W': 186.21, 'Y': 163.18, 'V': 99.13,
    }
    mass = sum(amino_acid_mass_dict[aa] for aa in peptide_seq)
    return {peptide_seq: mass}

def calculate_mol_mass_collection(peptides, amino_acid_mass_dict):
    """
     Accepts a list of peptides and returns a dictionary mapping each peptide to its neutral mass.
    
    Parameters
    ----------
    peptides : list
        A list of peptide sequences (strings).
    amino_acid_mass_dict : dict
        A dictionary mapping each amino acid to its corresponding mass in Daltons.

    Returns
    -------
    dict
        A dictionary mapping each peptide to its neutral mass.

    """
    peptide_masses = {}

    for peptide in peptides:
        peptide_masses[peptide] = calculate_mol_mass(
            peptide, amino_acid_mass_dict
        )[peptide]

    return peptide_masses

def calculate_mz_collection(peptide_mass_map, charge, proton_mass=1.007):
    """
    Accepts a peptide-mass dictionary and returns a new dictionary mapping peptides to their m/z values.

    Parameters
    ----------
    peptide_mass_map : dict
        A dictionary mapping each peptide to its neutral mass.
    charge : int
        The charge state to be used for m/z calculation.
    proton_mass : float, optional
        The mass of a proton in Daltons (default is 1.007).

    Returns
    -------
    dict
        A dictionary mapping each peptide to its m/z value.

    """
    mz_map = {}

    for peptide, neutral_mass in peptide_mass_map.items():
        mz_map[peptide] = (
            neutral_mass + charge * proton_mass
        ) / charge

    return mz_map

import numpy as np
import matplotlib.pyplot as plt

def plot_spectrum(mz_values, random_count_range=(0, 30000), seed=42):
    """
    Accepts the m/z values and optional `random_count_range` to draw random intensities, then plots a bar chart using Matplotlib.

    Parameters
    ----------
    mz_values : dict
        A dictionary mapping each peptide to its m/z value.
    random_count_range : tuple, optional
        The range of random counts to draw (default is (0, 30000)).
    seed : int, optional
        The seed for the random number generator (default is 42).

    Returns
    -------
    None
        This function plots the spectrum but does not return any value.

    -------

    """
    np.random.seed(seed)
    intensities = np.random.randint(
        random_count_range[0], random_count_range[1], size=len(mz_values)
    )

    plt.figure(figsize=(10, 6))
    plt.bar(mz_values.values(), intensities, width=1, color='blue', alpha=0.7)
    plt.xlabel('m/z')
    plt.ylabel('Intensity')
    plt.title('Simulated Mass Spectrum')
    plt.grid(axis='y')
    plt.show()

def fragment_peptide(peptide):
    """
    returns the combined list of b- and y-ion sequences for any peptide input.

    Parameters
    ----------
    peptide : str
        The amino acid sequence of the peptide to fragment.

    Returns
    -------
    list
        A list of the b- and y-ion sequences.

    """
    b_ions = [peptide[:i] for i in range(1, len(peptide))]
    y_ions = [peptide[i:] for i in range(1, len(peptide))]
    return b_ions + y_ions