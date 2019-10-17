""" A molecule """

import copy
import numpy

#import atom
from .atom import Atom
from .bond import Bond
from .angle import Angle

__has_openbabel__ = False
try:
    from openbabel import openbabel
except ImportError:
    pass
else:
    __has_openbabel__ = True


class BaseMolecule(object):
    """ A molecule

        A molecule is a collection of atoms and provides methods to
        extract or manipulate atoms.

        Internally, all atomic coordinates are stored in Angstrom
        and this convention should be adhered to when attempting to
        update the coordinates.

        Molecules from other APIs derive from this base class and
        implements all methods defined by this abstract class.
        See for example both Molecule and OBMolecule classes below.

    """
    _bond_threshold = 0.45 # Added threshold for bonds. Replicates openbabel

    def __init__(self):
        """ Initializes an empty molecule """
        self._charge = 0
        self._multiplicity = 1
        self._atoms = []
        self._bonds = []
        self._name = ""

    #
    # class methods
    #
    @classmethod
    def from_molecule(cls, m):
        """ Loads a molecule from another molecule by copying all data

            :param cls: class
            :type cls: BaseMolecule
            :param m: molecule
            :type m: BaseMolecule
        """
        M = cls()
        M.set_charge(m.get_charge())
        M.set_multiplicity(m.get_multiplicity())
        M.set_name(m.get_name())
        if m.get_num_atoms() > 0:
            M.add_atoms(*m.get_atoms())
        if m.get_num_atoms() > 1:
            M.add_bonds(*m.get_bonds())

        return M

    #
    # built-in modifications
    #
    def __len__(self):
        """ Returns the number of atoms in the molecule

            :rtype: int
        """
        return self.get_num_atoms()

    #
    # methods to add information to the molecule
    #
    def add_atom(self, _atom):
        """ Adds an atom to the molecule

            :param _atom: the atom to add to the molecule
            :type _atom: atom.Atom
        """
        raise NotImplementedError # pragma: no cover


    def add_atoms(self, *args):
        """ Adds multiple atoms to the molecule """
        for _atom in args:
            self.add_atom(_atom)


    def add_bond(self, _bond):
        """ Adds a bond to the molecule """
        raise NotImplementedError # pragma: no cover


    def add_bonds(self, *args):
        """ Adds multiple bonds to the molecule """
        for _bond in args:
            self.add_bond(_bond)

    #
    # getters and setters for various properties
    #
    def get_num_atoms(self):
        """ Returns the number of atoms in the molecule """
        raise NotImplementedError # pragma: no cover


    def get_atom(self, idx):
        """ Returns an atom based on its index

            Arguments:
            idx -- the index (from 0 to getNumAtoms() -1)

            Returns:
            an atom with the specified index
        """
        raise NotImplementedError # pragma: no cover


    def get_atoms(self):
        """ Returns an iterator of all atoms in the molecule """
        raise NotImplementedError # pragma: no cover


    def get_bonds(self):
        """ Returns an iterator of all bonds in the molecule """
        raise NotImplementedError # pragma: no cover


    def get_angles(self):
        """ Returns an iterator of all angles in the molecule """
        raise NotImplementedError # pragma: no cover

    #
    # getters and setters for simple properties
    #
    def get_name(self):
        """ Returns the name of the molecule """
        return self._name


    def set_name(self, value):
        """ Sets the name of the molecule

            Arguments:
            value -- the name of the molecule
        """
        assert isinstance(value, str)
        self._name = value


    def get_charge(self):
        """ Returns the integer charge of the molecule """
        return self._charge


    def set_charge(self, value):
        """ Sets the integer charge of the molecule

            Arguments:
            value -- the integer charge of the molecule
        """
        if not isinstance(value, int):
            raise ValueError("Argument 'value' must be of type integer")
        self._charge = value


    def get_multiplicity(self):
        """ Returns the multiplicity of the molecule """
        return self._multiplicity


    def set_multiplicity(self, value):
        """ Sets the multiplicity of the molecule

            Arguments:
            value -- the integer charge of the molecule
        """
        if not isinstance(value, int):
            raise ValueError("Argument 'value' must be of type integer")
        self._multiplicity = value

    #
    # properties that are lazily evaluated such as bonds
    #
    def percieve_bonds(self):
        """ Attempts to percieve covalent bonds in the molecule """
        raise NotImplementedError # pragma: no cover


    #
    # specialized getters and setters to extract information stored
    # in other classes related to the molecule such as atoms
    #
    def get_coordinates(self):
        """ Returns a numpy array with all the coordinates of the molecule

            Note: coordinates are always in Angstrom.
        """
        c = numpy.zeros((self.get_num_atoms(), 3))
        for iat, _atom in enumerate(self.get_atoms()):
            c[iat] = _atom.get_coordinate()

        return c


    def set_coordinates(self, value):
        """ Sets coordinates of molecule

            Arguments:
            value -- coordinates to store. Must be numpy array.
        """
        raise NotImplementedError # pragma: no cover


    def get_center_of_mass(self):
        """ Calculates the center of mass in units of Angstrom """
        mass = 0.0
        Rcm = numpy.zeros(3)
        for _atom in self.get_atoms():
            atom_mass = _atom.get_mass()
            mass += atom_mass
            Rcm += _atom.get_coordinate() * atom_mass

        assert mass != 0.0, "Total mass of molecule cannot be zero."

        return Rcm / mass

    #
    # Specialized iterators
    #
    def iter_atom_atoms(self, atom):
        """ Iterates over all atoms covalently bound to an atom

            Arguments:
            atom -- the atom whose neighbours to get
        """

        neighbour_indices = []
        for bond in self.get_bonds():
            try:
                nbr_index = bond.get_nbr_atom_idx(atom.get_idx())
            except ValueError:
                continue
            else:
                neighbour_indices.append(nbr_index)

        for _atom in self.get_atoms():
            if _atom.get_idx() in neighbour_indices:
                yield _atom


    def iter_atom_angles(self, atom):
        """ Iterates over all angles where atom is the vertex

            Arguments:
            atom -- the atom to use as the vertex in all angles
        """
        for angle in self.get_angles():
            if angle[0] == atom.get_idx():
                yield angle

    def iter_bond_atoms(self, bond):
        """ Iterates over all atoms in a bond

            Arguments:
            bond -- the bond to iterate over
        """
        for _atom in self.get_atoms():
            if _atom.get_idx() == bond._id1 or _atom.get_idx() == bond._id2:
                yield _atom


    def find_children(self, other_atom):
        """ Finds all atoms in the molecular graph as the supplied atom

            :param other_atom: the atom whose neighbours to get
            :type other_atom: atom.Atom
            :returns: a list of atoms sorted according to their internal index
            :rtype: list
        """
        raise NotImplementedError # pragma: no cover


class Molecule(BaseMolecule):
    """ A molecule

        A molecule is a collection of atoms and provides methods to
        extract or manipulate atoms.

        A library-free implementation of the Molecule class. The
        goal is provide a solution (albeit it might be slow) that
        does not depend on any third-party libraries.

        Internally, all atomic coordinates are stored in Angstrom
        and this convention should be adhered to when attempting to
        update the coordinates.
    """
    def __init__(self):
        """ Initializes an empty molecule """
        BaseMolecule.__init__(self)


    def add_atom(self, _atom):
        """ Adds an atom to the molecule """
        if not isinstance(_atom, Atom):
            raise TypeError
        self._atoms.append(copy.deepcopy(_atom))


    def add_bond(self, _bond):
        """ Adds a bond to the molecule

            :param _bond: The bond.Bond to add
            :type _bond: bond.Bond

            Checks are performed upon addition of bonds that it does not
            exist beforehand
        """
        #print("{}.{}({})".format(type(self).__name__, "add_bond", _bond))
        if not _bond in self._bonds:
            self._bonds.append(_bond)


    def get_num_atoms(self):
        """ Returns the number of atoms in the molecule

           :rtype: int
        """
        return len(self._atoms)


    def get_atom(self, idx):
        """ Returns an atom based on its index

            :param idx: the index (from 0 to get_num_atoms() -1)
            :type idx: int
            :rtype: atom.Atom
        """
        n_atoms = self.get_num_atoms()
        if idx < 0:
            raise IndexError("argument idx to getAtom must be >= 0")
        if idx >= n_atoms:
            raise IndexError("argument idx to getAtom must be < {}".format(n_atoms))
        return self._atoms[idx]


    def get_atoms(self):
        """ Returns all atoms (as an iterator) in the molecule

            :returns: all atoms as an iterator
            :rtype: collections.Iterable[atom.Atom]
        """
        for _atom in self._atoms:
            yield _atom


    def get_bonds(self):
        """ Returns an iterator of all bonds in the molecule

            If the bond list has not been calculated before, the bonds are
            percieved through the percieveBonds method

            :returns: all bonds as an iterator
            :rtype: collections.Iterable[bond.Bond]
        """
        if len(self._bonds) == 0:
            self._bonds = list(self.percieve_bonds())

        for _bond in self._bonds:
            yield _bond


    def get_angles(self):
        """ Returns an iterator of all angles in the molecule

            :returns: all angles as an iterator
            :rtype: collections.Iterable[angle.Angle]
        """
        for ibd, bond1 in enumerate(self.get_bonds()):
            for jbd, bond2 in enumerate(self.get_bonds()):
                if ibd <= jbd: continue
                jatm = bond1.shares_atom(bond2)
                if jatm >= 0:
                    iatm = bond1.get_nbr_atom_idx(jatm)
                    katm = bond2.get_nbr_atom_idx(jatm)
                    yield Angle(jatm, iatm, katm)


    def percieve_bonds(self):
        """ Attempts to percieve covalent bonds in the molecule

            It compares atom distances to covalent radii of the atoms.
        """
        #print("{}.{}".format(type(self).__name__, "percieve_bonds"))
        for iat, atom1 in enumerate(self.get_atoms()):
            for jat, atom2 in enumerate(self.get_atoms()):
                if iat <= jat: continue
                dr = atom2.get_coordinate() - atom1.get_coordinate()
                R2 = dr.dot(dr)

                dr_cov = atom1.get_covalent_radius() + atom2.get_covalent_radius() + self._bond_threshold
                R2_cov = dr_cov**2
                if R2 < R2_cov:
                    yield Bond(id1=atom1.get_idx(), id2=atom2.get_idx())


    def set_coordinates(self, value):
        """ Sets coordinates of molecule

            Arguments:
            value -- coordinates to store. Must be numpy array.
        """
        if not isinstance(value, numpy.ndarray):
            raise TypeError("Argument 'value' must be of type numpy array")
        (n, ) = numpy.shape(value)
        if n != self.get_num_atoms():
            raise ValueError("Argument 'value' has the wrong number of atoms")
        for iat, _atom in enumerate(self.get_atoms()):
            _atom.set_coordinate(value[iat])


    def find_children(self, other_atom):
        """ Finds all atoms in the molecular graph as the supplied atom

            :param other_atom: the atom whose neighbours to get
            :type other_atom: atom.Atom
            :returns: a list of atoms sorted according to their internal index
            :rtype: list
        """
        atoms = []
        old_atoms = [other_atom]

        for i in range(2):
            atoms = []
            for _atom in old_atoms:
                new_atoms = []
                for _nbr in self.iter_atom_atoms(_atom):
                    if _nbr not in old_atoms:
                        new_atoms.append(_nbr)

                if new_atoms: # check if sequency is empty
                    atoms.extend(new_atoms)

            if atoms: # check if sequency is empty
                old_atoms.extend(atoms)
            else:
                break

        return sorted(old_atoms, key=lambda _atom: _atom.get_idx())

class OBMolecule(BaseMolecule):
    def __init__(self, fromOBMol=None):
        raise NameError("Could not find OpenBabel.")

if __has_openbabel__:

    class OBMolecule(BaseMolecule):
        """ A molecule

            A molecule is a collection of atoms and provides methods to
            extract or manipulate atoms.

            This is an implementation of the Molecule class that directly
            supports openbabel as the underlying engine.

            Internally, all atomic coordinates are stored in Angstrom
            and this convention should be adhered to when attempting to
            update the coordinates.
        """
        def __init__(self, fromOBMol=None):
            """ Initializes an empty molecule """
            BaseMolecule.__init__(self)
            if fromOBMol is not None:
                self._obmol = openbabel.OBMol(fromOBMol)
                self.set_name(fromOBMol.GetTitle())
            else:
                self._obmol = openbabel.OBMol()


        def add_atom(self, _atom):
            """ Adds an atom to the molecule

                This method creates OBAtoms based on the atoms so it can use
                the internal functions that openbabel provides through OBMol.

                Note: This will surely break at some point for .pdb files etc.
            """
            if not isinstance(_atom, Atom):
                raise TypeError
            _obatom = openbabel.OBAtom()
            _obatom.SetAtomicNum(_atom.get_nuclear_charge())
            x, y, z = _atom.get_coordinate()
            _obatom.SetVector(x, y, z)
            _obatom.SetId(_atom.get_idx())
            _obatom.SetFormalCharge(_atom.get_formal_charge())
            self._obmol.AddAtom(_obatom)


        def add_bond(self, _bond):
            #print("{}.{}({})".format(type(self).__name__, "add_bond", _bond))
            _obatoms = []
            for _obatom in openbabel.OBMolAtomIter(self._obmol):
                #_atom = atom.Atom.from_obatom(_obatom)
                try:
                    ob_atom_index = int(_obatom.GetId()) # openbabel returns a long here for python2
                    atom_index = _bond.get_nbr_atom_idx(ob_atom_index)
                except ValueError:
                    pass
                else:
                    _obatoms.append(_obatom)

            assert len(_obatoms) == 2, "number of atoms in a bond must be 2."

            _obbond = openbabel.OBBond()
            _obbond.SetBegin(_obatoms[0])
            _obbond.SetEnd(_obatoms[1])
            _obbond.SetBondOrder(_bond.get_bond_order())
            self._obmol.AddBond(_obbond)


        def get_multiplicity(self):
            """ Returns the multiplicity of the molecule """
            return self._obmol.GetTotalSpinMultiplicity()


        def set_multiplicity(self, value):
            """ Sets the multiplicity of the molecule

                Arguments:
                value -- the integer charge of the molecule
            """
            if not isinstance(value, int):
                raise ValueError("Argument 'value' must be of type integer")
            self._obmol.SetTotalSpinMultiplicity(value)

        def get_charge(self):
            """ Returns the integer charge of the molecule """
            return self._obmol.GetTotalCharge()


        def set_charge(self, value):
            """ Sets the integer charge of the molecule

                Arguments:
                value -- the integer charge of the molecule
            """
            if not isinstance(value, int):
                raise ValueError("Argument 'value' must be of type integer")
            self._obmol.SetTotalCharge(value)


        def get_num_atoms(self):
            """ Returns the number of atoms in the molecule """
            return self._obmol.NumAtoms()


        def get_atom(self, idx):
            """ Returns an atom based on its index

                Note: The internal format for indices in openbabel is
                      to count from 1 to getNumAtoms() which is offset
                      by one compared to the format we have chosen in
                      the molecool library

                Arguments:
                idx -- the index (from 0 to getNumAtoms() -1)

                Returns:
                an atom with the specified index

            """
            n_atoms = self.get_num_atoms()
            if idx < 0:
                raise IndexError("argument idx to getAtom must be >= 0")
            if idx >= n_atoms:
                raise IndexError("argument idx to getAtom must be < {}".format(n_atoms))
            return Atom.from_obatom(self._obmol.GetAtom(idx+1))


        def get_atoms(self):
            """ Returns all atoms (as an iterator) in the molecule """
            for _obatom in openbabel.OBMolAtomIter(self._obmol):
                yield Atom.from_obatom(_obatom)


        def get_bonds(self):
            """ Returns an iterator of all bonds in the molecule """
            #self._obmol.ConnectTheDots() # what happens if called more than once?
            for _obbond in openbabel.OBMolBondIter(self._obmol):
                iat = _obbond.GetBeginAtom()
                jat = _obbond.GetEndAtom()
                yield Bond(iat.GetId(), jat.GetId(), order=_obbond.GetBondOrder())


        def get_angles(self):
            """ Returns an iterator of all angles in the molecule """
            self._obmol.FindAngles() # does nothing if angles have been found
            for _obangle in openbabel.OBMolAngleIter(self._obmol):
                vertex, id1, id2 = _obangle
                yield Angle(vertex, id1, id2)


        def percieve_bonds(self):
            """ Attempts to percieve covalent bonds in the molecule """
            self._obmol.ConnectTheDots()


        def set_coordinates(self, value):
            """ Sets coordinates of molecule

                Arguments:
                value -- coordinates to store. Must be numpy array.
            """
            if not isinstance(value, numpy.ndarray):
                raise TypeError("Argument 'value' must be of type numpy array")
            (n, ) = numpy.shape(value)
            if n != self.get_num_atoms():
                raise ValueError("Argument 'value' has the wrong number of atoms")
            for iat, _obatom in enumerate(openbabel.OBMolAtomIter(self._obmol)):
                x, y, z = value[iat]
                _obatom.SetVector(x, y, z)

        def find_children(self, other_atom):
            """ Finds all atoms in the molecular graph as the supplied atom

                :param other_atom: the atom whose neighbours to get
                :type other_atom: atom.Atom
                :returns: a list of atoms sorted according to their internal index
                :rtype: list
            """
            self._obmol.ConnectTheDots() # what happens if called more than once?
            idx = int(other_atom.get_idx()+1)
            vector = openbabel.vectorInt()
            self._obmol.FindChildren(vector, 0, idx)
            indices = sorted([i-1 for i in vector] + [idx-1])
            atoms = [self.get_atom(i) for i in indices]
            return sorted(atoms, key=lambda _atom: _atom.get_idx())
