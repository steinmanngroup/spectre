import spectre.errors

def get_chromophore_peex_data(filename, coupling_with_moments):
    """ Parses chromophore data from a DALTON log file

        Arguments:
        filename -- filename to parse from

        Returns:
        excitation energies
    """

    energies = []
    transition_dipoles = []
    tr_charges = []

    tr_dip = None

    parsing_eex_data = False

    with open("{0:s}.out".format(filename), "r") as peex_file:
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
                    tr_dip.append(float(tokens[9]))
                    if len(tr_dip) == 3:
                        transition_dipoles.append(tr_dip)

                    if len(tr_dip) > 3:
                        raise ValueError("wrong data parsed for transition dipole moments.")


            # reading of charges and multipole moments are separate from the above
            if coupling_with_moments and parsing_eex_data and "Potential fitted multipole moments (QFITLIB)" in line:
                charges = []

                for k in range(3):
                    line = peex_file.readline()

                # check this line has the correct information
                if not "Charges" in line:
                    print(line)
                    raise ValueError("Format of DALTON log file wrong. Expected 'Charges:' but got '{0:s}'".format(line))
                else:
                    line = peex_file.readline()

                # now parse charge data
                line = peex_file.readline()
                while "@" in line:
                    tokens = line.split()
                    charges.append(float(tokens[2]))
                    line = peex_file.readline()

                print(charges)

                tr_charges.append(charges)

            line = peex_file.readline()

    if len(transition_dipoles) == 0 and not coupling_with_moments:
        raise spectre.errors.SpectreExcitedStateValueError("No transition dipoles found in file '{0:s}'.".format(filename))

    if len(tr_charges) == 0 and coupling_with_moments:
        raise spectre.errors.SpectreExcitedStateValueError("No transition charges found in file '{0:s}'.".format(filename))

    return energies, transition_dipoles, tr_charges