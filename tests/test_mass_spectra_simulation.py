def test_fragment_peptide():
    peptide = "PEPTIDE"
    expected = ['P', 'PE', 'PEP', 'PEPT', 'PEPTI', 'PEPTID', 'EPTIDE', 'PTIDE', 'TIDE', 'IDE', 'DE', 'E']  # Example expected output
    actual = fragment_peptide(peptide)

    assert set(actual) == set(expected) 
    assert "PEPTIDE" not in expected  # which string should not be in the expected?

test_fragment_peptide()