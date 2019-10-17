import numpy

from spectre.errors import SpectreExcitedStateValueError

class SpectreExcitedStateData(object):
    """ Representation of excited state data in SPECTRE """

    def __init__(self):
        self._excitation_energies = []  # excitation energies
        self._tr_dips = []  # transition dipoles
        self.tr_q = []  # transition density fitted charges
        self.tr_d = []  # transition density fitted dipoles
        self.tr_o = []  # transition density fitted quadrupoles

    @classmethod
    def from_data(cls, de, tr_d, tr_moments, mom_order):
        """ Initiates the SpectreExcitedStateData class from data

            :param de: excitation energies
            :type de: list[float]
            :param tr_d: transition dipoles
            :type tr_d: list[list[float]]
            :param tr_moments: transition density fitted moments
            :type tr_moments: dict
            :param mom_order: the order of the multipole moments read from disk (0 = charges, 1 = charges, dipoles ...)
            :type mom_order: int

            :return: A populated SpectreExcitedStateData class
        """
        a = cls()
        #print("Created SpectreExcitedStateDate with", tr_moments)
        a.set_excitation_energies(de)
        a.set_transition_dipoles(tr_d)

        #print("mom_order:", mom_order)
        if mom_order not in [0, 1, 2]:
            raise ValueError("No multipole data included.")

        tr_q = tr_moments["charges"]
        nex, nat = numpy.shape(tr_q)
        a.set_transition_density_fitted_charges(tr_q)

        # we fill the arrays with blank data
        if mom_order >= 0:
            tr_d = numpy.zeros((nex, nat, 3))
            a.set_transition_density_fitted_dipoles(tr_d)

            tr_o = numpy.zeros((nex, nat, 6))
            a.set_transition_density_fitted_quadrupoles(tr_o)
        if mom_order >= 1:
            tr_d = tr_moments["dipoles"]
            a.set_transition_density_fitted_dipoles(tr_d)

            tr_o = numpy.zeros((nex, nat, 6))
            a.set_transition_density_fitted_quadrupoles(tr_o)

        if mom_order >= 2:
            tr_o = tr_moments["quadrupoles"]
            a.set_transition_density_fitted_quadrupoles(tr_o)
        return a

    def get_number_of_excited_states(self):
        return len(self._excitation_energies)

    def get_excitation_energies(self):
        if len(self._excitation_energies) == 0:
            raise SpectreExcitedStateValueError("No excitation energies stored.")

        return self._excitation_energies

    def set_excitation_energies(self, value):
        self._excitation_energies = numpy.array(value)

    def get_transition_dipoles(self):
        return self._tr_dips

    def set_transition_dipoles(self, value):
        self._tr_dips = numpy.array(value)

    def get_transition_density_fitted_charges(self):
        return self.tr_q

    def set_transition_density_fitted_charges(self, value):
        self.tr_q = numpy.array(value)

    def get_transition_density_fitted_dipoles(self):
        return self.tr_d

    def set_transition_density_fitted_dipoles(self, value):
        self.tr_d = numpy.array(value)

    def get_transition_density_fitted_quadrupoles(self):
        return self.tr_o

    def set_transition_density_fitted_quadrupoles(self, value):
        self.tr_o = numpy.array(value)

    def get_oscillator_strengths(self):
        """ Computes and returns the oscillator strengths of excitation

            The oscillator strenth is computed from

            :math:`2/3 \\Delta E |\\mu_T|^2`

            where :math:`\\Delta E` is an excitation energy and :math:`\\mu_T` is the
            transition dipole moment associated with this excitation energy.

            :return: oscillator strengths for all excitations
            :rtype: list[float]
        """
        osc = []
        for e, d in zip(self.get_excitation_energies(), self.get_transition_dipoles()):
            osc.append(2.0/3.0 * e * d.dot(d))

        return numpy.array(osc)

    def __str__(self):
        line = "@  1{0:>7d}{1:11.4f}{2:14.4f}{3:9.4f}{4:12.3f}{5:9.3f}{6:9.3f}\n"
        s = ""
        for i, (e, d) in enumerate(zip(self._excitation_energies, self._tr_dips), start=1):
            s += line.format(i, e, 0.0, 2./3.*e*d.dot(d), 0.0, 0.0, 0.0)
        return s