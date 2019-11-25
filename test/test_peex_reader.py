import spectre.readers


def test_get_chromophore_peex_data_no_charges():
    filename = "test/peex"
    e, dip, q, _ = spectre.readers.get_chromophore_peex_data(filename, False)
    assert len(e) == len(dip)


def test_get_chromophore_peex_data_with_charges():
    filename = "test/peex"
    e, dip, moments, _ = spectre.readers.get_chromophore_peex_data(filename, True)
    assert len(e) == len(moments['charges'])


if __name__ == '__main__':
    test_get_chromophore_peex_data_with_charges()
