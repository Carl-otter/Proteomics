def read_fasta(filepath):
    """
    Read a FASTA file and return a dictionary with Protein IDs attached to the sequence 

    Parameters
    ----------
    filepath : str
        Path to the FASTA file
    Returns
    -------
    dict
        Dictionary with the Protein IDs as keys for the corresponding sequence
    """

    fasta_path = filepath # Path to FASTA file as variablw
    protein_map = {}                            # Dictionary that will hold the protein IDs and sequences
    current_id = None                           # Variable that wl hold the current protein ID inside the loop
    current_sequence = []                       # List that will hold the current protein sequence inside the loop

    with open(fasta_path, 'r') as fasta_handle:                             # Opens the FASTA file for reading
        for line in fasta_handle:                                           # loop that goes over each line of the FASTA file 
            stripped = line.strip()                                         # Strips whitespace from the line
            if not stripped:                                                # If the line is fully empty it skips to the next line of the loop
                continue
            if stripped.startswith('>'):                                    # If the line starts with a ">" it is identified as a header line 
                if current_id is not None:                                  # If the current_id is not None, this means that we have already read a protein sequence and we can add it to the dictionary
                    protein_map[current_id] = ''.join(current_sequence)     # Joins the protein sequence lines to a single line and adds them to the dictionary with the ID as key
                    current_sequence = []                                   # Resets the current_sequence list for the next protein sequence  
                current_id = stripped.split('|')[1]                         # Extracts the protein ID from the header line and assigns it to current_id
            else:
                current_sequence.append(stripped)                           # if the line is not a headerline the line is added to the current sequence list

    if current_id is not None:                                              # if the current id is not none, this means that we have read a protein sequence and we can add it to the dictionary (end of loop)
        protein_map[current_id] = ''.join(current_sequence)                 # joins the protein sequence to a single line and adds them to the dictionary

    return protein_map                                                             # returns the dictionary

    