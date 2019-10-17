import numpy
import os
import os.path

from spectre.errors import SpectrePEEXFileNotFoundError, SpectreExcitedStateValueError


def get_chromophore_peex_data(filename, coupling_with_moments):
    """ Parses chromophore data from a DALTON log file

        :param filename: the file to read
        :type filename: str
        :param coupling_with_moments: whether or not to read the coupling data
        :type coupling_with_moments: bool
        :raises: if file is not found SpectrePEEXFileNotFoundError is twrown.
        :return: a lot of data
        :rtype: tuple[list[float], list[list[float]], dict, int]
    """

    energies = []
    transition_dipoles = []
    tr_moments = {}

    tr_charges = []
    tr_dipoles = []
    tr_quadrupoles = []
    if coupling_with_moments:
        tr_moments["charges"] = []
        tr_moments["dipoles"] = []
        tr_moments["quadrupoles"] = []

    tr_dip = None

    parsing_eex_data = False
    parsing_tr_data = False
    peex_filename = "{0:s}.out".format(filename)
    if not os.path.exists(peex_filename):
        raise SpectrePEEXFileNotFoundError("Could not find the excited state file {}".format(peex_filename))
    with open(peex_filename, "r") as peex_file:
        line = peex_file.readline()

        while line:

            # we start by reading an excitation energy from the DALTON log file
            if "@ Excitation energy :" in line:
                tokens = line.split()
                energies.append(float(tokens[4]))

                tr_dip = []  # storage for the upcoming transition dipole moments
                parsing_eex_data = True  # flag we are parsing data now

            if tr_dip is not None:
                if "@ Oscillator strength (LENGTH)" in line:
                    tokens = line.split()
                    tr_dip.append(float(tokens[5]))
                    if len(tr_dip) == 3:
                        transition_dipoles.append(tr_dip)

                    if len(tr_dip) > 3:
                        raise ValueError("wrong data parsed for transition dipole moments.")

            # reading of charges and multipole moments are separate from the above
            if coupling_with_moments and parsing_eex_data and "Potential fitted multipole moments (QFITLIB)" in line:
                parsing_tr_data = True

                # now we start to read moments
            if "Charges:" in line and "charges" in tr_moments and parsing_tr_data:
                charges = []
                peex_file.readline()  # the line with Q
                line = peex_file.readline()  # the first data line
                while "@" in line:
                    tokens = line.split()
                    charges.append(float(tokens[2]))
                    line = peex_file.readline()

                tr_charges.append(charges)

            if "Dipoles:" in line and "dipoles" in tr_moments and parsing_tr_data:
                dipoles = []
                peex_file.readline()  # the line with Q
                line = peex_file.readline()  # the first data line
                while "@" in line:
                    tokens = line.split()
                    dipoles.append(list(map(float, tokens[2:])))
                    line = peex_file.readline()

                tr_dipoles.append(dipoles)

            if "Quadrupoles:" in line and "quadrupoles" in tr_moments and parsing_tr_data:
                quads = []
                peex_file.readline()  # the line with Q
                line = peex_file.readline()  # the first data line
                while "@" in line:
                    tokens = line.split()
                    # expected format is
                    # XX          XY          XZ          YY          YZ          ZZ
                    quads.append(list(map(float, tokens[2:])))
                    line = peex_file.readline()

                tr_quadrupoles.append(quads)

            line = peex_file.readline()

    if len(transition_dipoles) == 0 and not coupling_with_moments:
        raise SpectreExcitedStateValueError("No transition dipoles found in file '{0:s}'.".format(filename))

    if coupling_with_moments:
        tr_moments["charges"] = numpy.array(tr_charges)
        tr_moments["dipoles"] = numpy.array(tr_dipoles)
        tr_moments["quadrupoles"] = numpy.array(tr_quadrupoles)
        if len(tr_moments["charges"]) == 0:
            raise SpectreExcitedStateValueError("No transition charges found in file '{0:s}.out'.".format(filename))

    mom_order = -1
    if coupling_with_moments:
        if len(tr_moments["charges"]) > 0:
            mom_order = 0
        if len(tr_moments["dipoles"]) > 0:
            mom_order = 1
        if len(tr_moments["quadrupoles"]) > 0:
            mom_order = 2
    return energies, transition_dipoles, tr_moments, mom_order


def test_m0_reader():
    filename = "m0"
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
    filename = "m1"
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
    filename = "m2"
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