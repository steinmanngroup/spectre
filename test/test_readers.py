from spectre.readers import get_chromophore_peex_data


def test_m0_reader():
    filename = "test/m0"
    erg1, trd1, trm1, mo = get_chromophore_peex_data(filename, False)
    assert type(erg1) == list
    assert len(erg1) == 4
    assert len(erg1) == len(trd1)
    assert type(trm1) == dict
    assert mo == -1  # no moments
    assert abs(trd1[3][2] - 2.08899512e-02) < 1.0e-9  # checking on the last item in the list of transition dipoles

    # now read again, but this time read the m0 moments too
    erg2, trd2, trm2, mo = get_chromophore_peex_data(filename, True)
    assert len(erg2) == len(erg1)
    for e1, e2 in zip(erg1, erg2):
        assert abs(e2 - e1) < 1.0e-9
    assert type(trm2) == dict
    assert mo == 0  # charges
    assert "charges" in trm2
    tr_q = trm2["charges"]
    assert len(tr_q) == len(erg2)
    assert abs(tr_q[2][0] - (-0.566597)) < 1.0e-9  # the 3rd excitation first charge


def test_m1_reader():
    filename = "test/m1"
    erg1, trd1, trm1, mo = get_chromophore_peex_data(filename, False)
    assert type(erg1) == list
    assert len(erg1) == 4
    assert len(erg1) == len(trd1)
    assert type(trm1) == dict
    assert mo == -1  # no moments
    assert abs(trd1[0][0] - 1.50505287e-02) < 1.0e-9  # checking on the first item in the list of transition dipoles

    # now read again, but this time read the m0 and m1 moments too
    erg2, trd2, trm2, mo = get_chromophore_peex_data(filename, True)
    assert mo == 1  # dipoles
    assert "charges" in trm2
    tr_q = trm2["charges"]
    assert len(tr_q) == len(erg2)
    assert abs(tr_q[0][1] - 0.003357) < 1.0e-9  # the 1st excitation, 2nd atom charge

    assert "dipoles" in trm2
    tr_d = trm2["dipoles"]
    assert len(tr_d) == len(erg2)
    for i in range(len(tr_d)):
        assert len(tr_d[i]) == 3  # number of atoms
        for j in range(len(tr_d[i])):
            assert len(tr_d[i][j]) == 3
    assert abs(tr_d[1][1][1] - 0.013144) < 1.0e-9  # 2nd excitation, 2nd atom, 2nd element of dipole


def test_m2_reader():
    filename = "test/m2"
    erg1, trd1, trm1, mo = get_chromophore_peex_data(filename, False)
    assert type(erg1) == list
    assert len(erg1) == 4
    assert len(erg1) == len(trd1)
    assert type(trm1) == dict
    assert mo == -1  # no moments
    assert abs(trd1[2][1] - 6.28123166e-02) < 1.0e-9  # checking on the first item in the list of transition dipoles

    # now read again, but this time read the m0 and m1 moments too
    erg2, trd2, trm2, mo = get_chromophore_peex_data(filename, True)
    assert mo == 2  # quadrupoles
    assert "charges" in trm2
    tr_q = trm2["charges"]
    assert len(tr_q) == len(erg2)
    assert abs(tr_q[0][1] - (-0.473546)) < 1.0e-9  # the 1st excitation, 2nd atom charge

    assert "dipoles" in trm2
    tr_d = trm2["dipoles"]
    assert len(tr_d) == len(erg2)
    assert abs(tr_d[1][2][2] - 0.163796) < 1.0e-9  # 2nd excitation, 3nd atom, 3nd element of dipole

    assert "quadrupoles" in trm2
    tr_quad = trm2["quadrupoles"]
    assert len(tr_quad) == len(erg2)
    for i in range(len(tr_quad)):
        assert len(tr_quad[i]) == 3  # number of atoms
        for j in range(len(tr_quad[i])):
            assert len(tr_quad[i][j]) == 6
    assert abs(tr_quad[2][0][5] - (-0.279883)) < 1.0e-9  # 3rd excitation, 1st atom, 6th element of quad tensor
    assert abs(tr_quad[3][0][0] - 0.266974) < 1.0e-9  # 3rd excitation, 1st atom, 6th element of quad tensor
