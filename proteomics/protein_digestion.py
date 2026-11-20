def digest_protein_sequence(protein_seq, cleave_pattern, min_pep_len=5, max_pep_len=30):
    """
    The Preotein sequences are filtered to mimic experimental biased and then splitted at the celeavage site of the specific enzyme. 
    The function returns a list of the resulting peptide fragments.

    Parameters
    ----------
    protein_seq : str
        The amino acid sequence of the protein to be digested
    cleave_pattern : str
        The regex pattern that defines the cleavage site for the specific enzyme
    min_pep_len : int, optional
        The minimum length of the resulting peptide fragments (default is 5)
    max_pep_len : int, optional
        The maximum length of the resulting peptide fragments (default is 30)

    Returns
    -------
    list
        List of peptide fragments resulting from the digestion of the protei sequence
    """
    import re

    peptides = re.split(cleave_pattern, protein_seq)

    filtered = [peptide for peptide in peptides if min_pep_len <= len(peptide) <= max_pep_len]

    return filtered

def digest_protein_collection(protein_map, cleave_pattern, min_pep_len=5, max_pep_len=30):
    """
    Digests a cllection of proteins using a specific enzyme

    Parameters
    ----------
    protein_map : dict
        Dictionary with the Protein IDs as keys for the correspending sequence
    cleavage_pattern : str
        The regex pattern that defines the cleavage site for the specific enzyme
    min_pep_len : int, optional
        The minimum length of the resulting peptide fragments (default is 5)
    max_pep_len : int, optional
        The maximum length of the resulting peptide fragments (default is 30)

    Returns
    -------
    dict
        A dictionary containing the petide fragments of all digested proteins, with the Protein IDs as keys and the corresponding list of peptide fragments as values
    """
    digest_results = {}
    for protein_id, sequence in protein_map.items():
        peptides = digest_protein_sequence(sequence, cleave_pattern, min_pep_len, max_pep_len)
        digest_results[protein_id] = peptides
    return digest_results

def compute_sequence_coverage(protein_seq, peptides):
    """
    Return the percentage of protein residues covered by the peptides. Test if the peptide is located in the sequence.

    Parameters
    ----------
    protein_seq : str
        The amino acid sequence of the protein
    peptides : list
        List of peptide fragments resulting from the digestion of the protein sequence

    Returns
    -------
    float
        The percentage of protein residues covered by the peptides
    """

    covered_positions = set()

    for peptide in peptides:
        start = 0

        # Find every occurrence of this peptide in the protein.
        while peptide:
            start = protein_seq.find(peptide, start)
            if start == -1:
                break

            covered_positions.update(range(start, start + len(peptide)))
            start += 1

    return len(covered_positions) / len(protein_seq) * 100