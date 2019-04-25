import spectre.readers

def test_get_chromophore_peex_data_no_charges():
    filename = "peex"
    e, dip, q = spectre.readers.get_chromophore_peex_data(filename, False)
    assert len(e) == len(dip)

def test_get_chromophore_peex_data_with_charges():
    filename = "peex"
    e, dip, q = spectre.readers.get_chromophore_peex_data(filename, True)
    assert len(e) == len(q)

if __name__ == '__main__':
    test_get_chromophore_peex_data_with_charges()
