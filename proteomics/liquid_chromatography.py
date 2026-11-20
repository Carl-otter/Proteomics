from pyteomics import achrom

def predict_lc_retention_times(peptides):
    """
    Takes a list of peptide sequences and returns a dictionary mapping each peptide to its predicted retention time (rounded to two decimals)

    Parameters
    ----------
    peptides : list of str
        A list of peptide sequences.

    Returns
    -------
    dict of str -> float
        A dictionary mapping each peptide to its predicted retention time.
    
    Dependencies
    ------------
    Requires the pyteomics package
    """
    retention_times = {}
    for peptide in peptides:
        rt = achrom.calculate_RT(peptide, achrom.RCs_guo_ph7_0)
        retention_times[peptide] = round(rt, 2)
    return retention_times

import matplotlib.pyplot as plt

def plot_retention_time(retention_times, resolution=30):
    """
    Plots a histogram of the predicted retention times.

    Parameters
    ----------
    retention_times : dict
        A dictionary mapping each peptide to its predicted retention time.
    resolution : int, optional
        The number of bins in the histogram. Default is 30.
    Returns
    -------
    None
    """
    plt.figure(figsize=(10, 6))
    plt.hist(retention_times.values(), bins=resolution, color='blue', alpha=0.7)
    plt.title('Predicted Retention Time Distribution')
    plt.xlabel('Retention Time (minutes)')
    plt.ylabel('Number of Peptides')
    plt.grid(axis='y', alpha=0.75)
    plt.show()

def select_retention_time_window(peptide_rt_map, lower_ret_time, upper_ret_time):
    """
    Accept the retention-time dictionary and lower/upper bounds (in minutes) and return only the peptides that fall inside that window.

    Parameters
    ----------
    peptide_rt_map : dict
        A dictionary mapping each peptide to its predicted retention time.
    lower_ret_time : float
        The lower bound of the retention time window.
    upper_ret_time : float
        The upper bound of the retention time window.
    Returns
    -------
    dict
        A dictionary containing only the peptides that fall within the specified retention time window.

    """
    selected_peptides = {peptide: rt for peptide, rt in peptide_rt_map.items() if lower_ret_time <= rt <= upper_ret_time}
    return selected_peptides