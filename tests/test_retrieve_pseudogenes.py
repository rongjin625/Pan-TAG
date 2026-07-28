from pantag.retrieve_pseudogenes import merge_intervals


def test_merge_intervals_uses_union_length():
    assert merge_intervals([(1, 10), (5, 15), (20, 25)]) == 21
    assert merge_intervals([]) == 0
