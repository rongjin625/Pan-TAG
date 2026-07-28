from pantag.single_sample_tag import family_prefix, filter_blastp


def test_family_prefix_matches_legacy_rule():
    assert family_prefix("PRAMEF12") == "PRAMEF"
    assert family_prefix("USP17L") == "USP"


def test_filter_blastp_uses_identity_and_qcovhsp(tmp_path):
    raw = tmp_path / "hits.raw"
    output = tmp_path / "hits.filtered"
    raw.write_text(
        "q1\ts1\t25\t10\t0\t0\t1\t10\t1\t10\t1e-9\t50\t30\n"
        "q2\ts2\t19\t10\t0\t0\t1\t10\t1\t10\t1e-9\t50\t30\n"
        "q3\ts3\t25\t10\t0\t0\t1\t10\t1\t10\t1e-9\t50\t19\n",
        encoding="utf-8",
    )
    assert filter_blastp(raw, output, 20, 20) == 1
    assert output.read_text(encoding="utf-8").startswith("q1\ts1\t25")
