from pantag.merge import clean_values, overlap, raw_values


def test_names_are_cleaned_but_physical_ids_are_preserved():
    assert clean_values("PRAMEF_2,USP17L(3),None") == {"PRAMEF", "USP17L"}
    assert raw_values("gene_2,None") == {"gene_2"}


def test_overlap_matches_pipeline8_directional_score():
    forward, reverse, score = overlap({"A", "B"}, {"A", "B", "C", "D"})
    assert forward == 1.0
    assert reverse == 0.5
    assert score == 1.0
