""" An atom in a molecule """

import numpy

from .util import MASSES, VDWRADII, COVALENTRADII, COORDINATION, Z2LABEL

class Atom(object):
    """ An atom

        The minimum amount of information required is the
        nuclear charge Z of the atom.

        A lot of additional default information is added to the atom
        based on the nuclear charge if not provided by the instantia-
        tor. The information is:

          * mass
          * Van der Waal radius
          * covalent radius
          * coordination number
          * label
          * hybridization (sp, sp2 or sp3)

        Arguments:
        Z -- nuclear charge of atom. This argument is mandatory.

        Keyword Arguments:
        mass -- the mass of the atom in atomic units. Default is specified by using the nuclear charge.
        xyz -- the Cartesian coordinate of the atom in Angstrom. Default is origo.
        idx -- the atom index. If not specified, a value of -1 is assigned.
        fcharge -- formal charge used for cat- and anions. Default 0.
    """
    def __init__(self, Z, **kwargs):
        assert Z > 0, "Nuclear charge of atom must be greater than zero."
        self._z = Z
        self._c = numpy.array(kwargs.get('xyz', [0, 0, 0]))
        self._mass = kwargs.get('mass', MASSES[Z])
        self._idx = kwargs.get('idx', -1)
        self._fcharge = kwargs.get('fcharge', 0)
        self._vdw_radius = kwargs.get('vwdradius', VDWRADII[Z])
        self._cov_radius = kwargs.get('covradius', COVALENTRADII[Z])
        self._coordination = kwargs.get('coordination', COORDINATION[Z])
        self._hybridization = kwargs.get('hybridization', 0)
        self._label = Z2LABEL[Z]


    @classmethod
    def from_atom(cls, other):
        """ Creates an atom from another Atom

            :param other: the other atom
            :type other: Atom
        """
        atom = cls(other.get_nuclear_charge(),
                   mass=other.get_mass(),
                   fcharge=other.get_formal_charge(),
                   xyz=other.get_coordinate(),
                   idx=other.get_idx(),
                   hybridization=other.get_hybridization())
        return atom


    @classmethod
    def from_obatom(cls, _obatom):
        """ Constructs an Atom from an openbabel OBAtom

            :param _obatom: an openbabel OBAtom
            :type _obatom: openbabel.OBAtom
        """
        x, y, z = _obatom.GetX(), _obatom.GetY(), _obatom.GetZ()
        _atom = cls(_obatom.GetAtomicNum(),
                    xyz=[x, y, z],
                    idx=_obatom.GetId(), # not using internal index (GetIdx) from openbabel because it gets remapped
                    fcharge=_obatom.GetFormalCharge(),
                    hybridization=_obatom.GetHyb())
        return _atom


    def get_mass(self):
        """ Returns the mass of the atom """
        return self._mass


    def get_nuclear_charge(self):
        """ Returns the nuclear charge of the atom """
        return self._z


    def get_formal_charge(self):
        """ Returns the formal charge of the atom """
        return self._fcharge


    def set_formal_charge(self, value):
        """ Sets the integer formal charge of the atom

            :param value: the integer formal charge
            :type value: int
        """
        if not isinstance(value, int):
            raise TypeError

        self._fcharge = value


    def get_idx(self):
        """ Returns the internal index of the atom """
        return self._idx


    def set_idx(self, value):
        """ Sets the internal index of the atom

            :param value: the integer index
            :type value: int
        """
        if not isinstance(value, int):
            raise TypeError
        self._idx = value


    def get_label(self):
        """ Returns the human readable atomic label """
        return self._label


    def get_coordinate(self):
        """ returns the 3D coordinate of the atom """
        return self._c


    def set_coordinate(self, value):
        """ Sets the 3D coordinate of the atom

            :param value: the coordinate given as a numpy array
            :type value: numpy.ndarray
        """
        if not isinstance(value, numpy.ndarray):
            raise TypeError("Argument 'value' must be of type numpy array")
        (n, ) = numpy.shape(value)
        if n != 3:
            raise ValueError("Dimensions of data do not match. Expected 3 but got {}".format(n))
        self._c = value


    def get_vdw_radius(self):
        """ Returns the Van der Waal radius of the atom """
        return self._vdw_radius


    def get_covalent_radius(self):
        """ Returns the covalent radius of the atom """
        return self._cov_radius


    def set_coordination(self, value):
        """ Sets the coordination number

            :param value: the coordination number
            :type value: int
        """
        if not isinstance(value, int):
            raise TypeError

        max_coordination = molecool.util.COORDINATION[self._z]
        if value > max_coordination:
            raise ValueError("Coordination number too large.")

        self._coordination = value


    def get_coordination(self):
        """ Returns the coordination number """
        return self._coordination


    def set_hybridization(self, value):
        """ Sets the hybridization of the atom

            typical values are:
               1 = sp
               2 = sp2
               3 = sp3

            Arguments:
            value -- the sp hybridization of the atom
        """
        if not isinstance(value, int):
            raise TypeError("Argument 'value' must be of type integer.")
        if value < 1:
            raise ValueError("Hybridization equal to or below zero not valid.")

        self._hybridization = value

    def get_hybridization(self):
        """ Returns the hybridization of the atom

            Typical values return are:
               1 = sp
               2 = sp2
               3 = sp3
        """
        return self._hybridization


    def __eq__(self, other):
        """ Tests that two atoms are equal using various measures

            We test the following:
                * do the atoms have the same nuclear charge?

            Other things that could be tested for in the future or
            through another method implying stronger equivalence:
                * hybridization

            We DO NOT test the following:
                * are they placed ontop of each other

            :param other: the other atom
            :type other: Atom
        """
        if self.get_nuclear_charge() != other.get_nuclear_charge():
            return False

        #EPS = 1.0e-6
        #dr = self.get_coordinate() - other.get_coordinate()
        #R2 = dr.dot(dr)
        #return numpy.sqrt(R2) < EPS
        return True


    def __repr__(self):
        fmt_str = "Atom({0:d}, xyz=[{1[0]:.7f}, {1[1]:.7f}, {1[2]:.7f}, idx={2:d}])"
        return fmt_str.format(self.get_nuclear_charge(),
                              self.get_coordinate(),
                              self.get_idx())
