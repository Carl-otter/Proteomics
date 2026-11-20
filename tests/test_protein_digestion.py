from proteomics import digest_protein_sequence

def test_digest_protein_sequence():

    dummy_seq = "AAAKBBBBKCCCCKDDDD"

    peptides = digest_protein_sequence(
        dummy_seq,
        cleave_pattern="(?<=K)",
        min_pep_len=1,
        max_pep_len=30
    )

    expected = ["AAAK", "BBBBK", "CCCCK", "DDDD"]

    assert peptides == expected

from proteomics.your_module import digest_protein_collection


def test_digest_protein_collection():

    dummy_proteins = {
        "protein1": "AAAKBBBBBKCCCCCCKDDDD",
        "protein2": "EEEEEKFFFFFFKGGGGGG"
    }

    digested_peptides_map = digest_protein_collection(
        dummy_proteins,
        cleave_pattern="(?<=K)",
        min_pep_len=5,
        max_pep_len=30
    )

    expected = {
        "protein1": ["BBBBBK", "CCCCCCK"],
        "protein2": ["EEEEEK", "FFFFFFK", "GGGGGG"]
    }

    assert digested_peptides_map == expected