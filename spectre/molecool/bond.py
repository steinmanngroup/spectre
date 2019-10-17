""" Bond between atoms """

class Bond(object):
    """ A bond between two Atom objects.

        Currently, there is no reference to the actual Atom objects but instead
        the atoms are indexed in the Molecule class storing these bond indices.

        Bonds between two atoms `A` and `B`, i.e. `A` -- `B`, are constructed
        manually as follows.

        >>> b = Bond(A.get_idx(), B.get_idx())

        Usually this is not needed as the Molecule class contains data
        to automatically construct bonds between atoms.
    """
    def __init__(self, id1, id2, order=1):
        """ Instantiates the Bond class

        :param id1: the first atom index
        :type id1: int
        :param id2: the second atom index
        :type id2: int
        :param order: the bond order. Default is 1 which means a single bond
        :type order: int

        :rtype: Bond
        """
        self._id1 = id1
        self._id2 = id2
        self._bond_order = order # 1 single, 2 double, 3 triple, 4 any, 5 aromatic # i guess?
        if self._id1 == -1:
            raise ValueError("Bond Error: first index must be different from -1")
        if self._id2 == -1:
            raise ValueError("Bond Error: second index must be different from -1")
        if self._id1 == self._id2:
            raise ValueError("Bond Error: indices cannot refer to same atom")


    def shares_atom(self, other):
        """ Returns an atom index if two bonds shares an atom.
            If an atom is not found or the bonds are the same a -1 is returned.

            :param other: The bond to check to see if they share atoms.
            :type other: Bond
            :return: integer of atom the bonds share. -1 if None.
            :rtype: int
        """
        if self == other:
            return -1

        if self._id1 == other._id1 and self._id2 != other._id2:
            return self._id1

        if self._id1 == other._id2 and self._id2 != other._id1:
            return self._id1

        if self._id2 == other._id1 and self._id1 != other._id2:
            return self._id2

        if self._id2 == other._id2 and self._id1 != other._id1:
            return self._id2

        return -1


    def get_bond_order(self):
        """ Returns the bond order of the bond

            :rtype: int
        """
        return self._bond_order


    def get_nbr_atom_idx(self, atom_index):
        """ Returns the neighboring atom index in the bond

            :param atom_index: The atom index to find the neighbour to
            :type atom_index: int
        """
        if not isinstance(atom_index, int):
            raise TypeError("Expected integer in get_nbr_atom_idx in class Bond. Got '{}' instead.".format(type(atom_index)))
        if self._id1 == atom_index:
            return self._id2
        if self._id2 == atom_index:
            return self._id1
        raise ValueError("The atom index {0:d} is not in the bond.".format(atom_index))


    def __eq__(self, other):
        """ Bonds are equal if they refer to the same atoms

            :param other: the other atom to compare to
            :type other: Bond
        """
        check_1 = self._id1 == other._id1 and self._id2 == other._id2
        check_2 = self._id1 == other._id2 and self._id2 == other._id1
        return check_1 or check_2


    def __repr__(self):
        return "Bond({0:d}, {1:d}, order={2:d})".format(self._id1, self._id2, self._bond_order)
