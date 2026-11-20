def test_read_fasta(tmp_path):
    from proteomics import read_fasta
    # Create a temporary FASTA file
    fasta_content = """>sp|Protein1|Some protein description
    GOOD
    BYE
    >sp|Protein2|Another protein description
    HEL
    LO
    """

    fasta_file = tmp_path / "test.fasta"
    fasta_file.write_text(fasta_content)

    # Run the function
    protein_map = read_fasta(fasta_file)

    # Check the expected output
    expected = {
        "Protein1": "GOODBYE",
        "Protein2": "HELLO"
    }

    assert protein_map == expected