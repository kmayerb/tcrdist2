import pytest
import os
from tcrdist import setup_db
tempSkip = pytest.mark.skip(reason="Temporarily skipping for efficiency.")
test_data = [('new_nextgen_chains_mouse_A.tsv', 'https://www.dropbox.com/s/pkpr6p97eworn3q/new_nextgen_chains_mouse_A.tsv?dl=1', None),
             ('new_nextgen_chains_mouse_B.tsv', 'https://www.dropbox.com/s/sxgvrj25mnzr20s/new_nextgen_chains_mouse_B.tsv?dl=1', None),
             ('new_nextgen_chains_human_A.tsv', 'https://www.dropbox.com/s/41w8yl38nr4ey32/new_nextgen_chains_human_A.tsv?dl=1', None),
             ('new_nextgen_chains_human_B.tsv', 'https://www.dropbox.com/s/8ysciqrcywdsryp/new_nextgen_chains_human_B.tsv?dl=1', None)
         ]
@tempSkip
@pytest.mark.parametrize("input, expect, expected_error", test_data)
def test_setub_db_nextgen_data_to_db(input, expect, expected_error):
    try:
        x = setup_db.install_nextgen_data_to_db(download_file = input)
        assert x.split(" ")[3] == expect
        assert os.path.isfile(x.split(" ")[2])
    except:
        with pytest.raises(expected_error):
            x = setup_db.install_nextgen_data_to_db(download_file = input)
