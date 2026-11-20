from proteomics.liquid_chromatography import predict_lc_retention_times
from pyteomics import achrom

def test_predict_lc_retention_times():

    peptides = ["PEPTIDE", "SEQUENCE", "PROTEIN"]

    retention_times = predict_lc_retention_times(peptides)

    expected = {
        peptide: round(
            achrom.calculate_RT(peptide, achrom.RCs_guo_ph7_0), 2
        )
        for peptide in peptides
    }

    assert retention_times == expected

from proteomics.protein_digestion import select_retention_time_window


def test_select_retention_time_window():

    peptide_rt_map = {
        "PEPTIDE": 5.5,
        "SEQUENCE": 12.3,
        "PROTEIN": 18.7,
        "AMINOACID": 25.4,
        "DIGESTION": 32.1
    }

    selected_peptides = select_retention_time_window(
        peptide_rt_map,
        lower_ret_time=10.0,
        upper_ret_time=25.0
    )

    expected = {
        "SEQUENCE": 12.3,
        "PROTEIN": 18.7
    }

    assert selected_peptides == expected