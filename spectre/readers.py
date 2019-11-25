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


