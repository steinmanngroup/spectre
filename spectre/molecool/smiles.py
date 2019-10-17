""" A SMILES engine """
from __future__ import print_function

from collections import Counter
import pytest

from molecool.molecule import Molecule
from molecool.atom import Atom
from molecool.bond import Bond
from molecool.util import LABEL2Z


BOND_ORDERS = {"-": 1, "=": 2, "#": 3, ":": 4}
MAX_RING_INDICES = 9


class AtomError(ValueError):
    pass


class AtomBracketError(AtomError):
    pass


class IllegalAtomError(AtomError):
    pass


class BranchBracketError(ValueError):
    pass


class RingError(ValueError):
    pass


class RingNumberNotSupported(RingError):
    pass


class RingClosureError(RingError):
    pass


def count_items_exclusive(s, values):
    """ counts occurrences of items in input values """
    counts = [0 for v in values]
    for i, c in enumerate(s):
        if c in values:
            ii = values.index(c)
            counts[ii] += 1

    return counts


class Smiles(object):
    """ The molecool SMILES engine

        Typical usage is::

        >>> S = Smiles("CC")

        which will generate a molecular graph, i.e. a molecule
        without the 3D geometry but with atoms and bonds like
        the following::

           Atom(6, idx=0) Bond(0, 1) Atom(6, idx=1)

        The engine supports the organic subset of atoms and
        atoms in brackets, i.e. [Fe2+] for example. There
        is support for rings (not nested though) and only
        up to 9 which means the % modifier has not been
        implemented yet.

        There is also support for branches.
    """
    def __init__(self, s, **kwargs):
        """ Initializes the SMILES engine

            :param s: the smiles string
            :type s: str
            :param kwargs: keyword arguments
            :type kwargs: dict

            :rtype: Smiles

            keyword arguments:
            debug -- enable debug printing. default is False.
        """
        self._mol = Molecule()

        self._debug = kwargs.get('debug', False)

        # vector of indices. A ring is present in atom i if
        # ring_openings[i] corresponds to an atom index in
        # the molecule
        self._ring_openings = [-1]*MAX_RING_INDICES

        if self._debug:
            print("------------------------------")
            print(s)
            print("------------------------------")
        self.parse_smiles(s, None)

        # charge of molecule = sum of formal charges on atoms as a first guess
        formal_charges = [a.get_formal_charge() for a in self._mol.get_atoms()]
        self._mol.set_charge(sum(formal_charges))

        # now we do some checks
        if sum(self._ring_openings) != -1*MAX_RING_INDICES:
            print(self._ring_openings)
            raise RingClosureError()

        if self._debug:
            print("------------------------------")


    def parse_smiles(self, s, p_atom, branch_index_offset=0):
        """ Parses (part of) a SMILES string

            Arguments:
            s -- the smiles string to parse
            p_atom -- the previous atom, if any, for bonding. Otherwise None.
            branch_index_offset -- an offset used when branching to update indices correctly.

            Returns:
            an updated atom index of the next atom
        """
        p_prev_atom = p_atom
        p_bond = None

        l_bracket_open = -1
        l_bracket_close = -1
        i_bracket_counter = 0 # ( and ) counters

        l_atom_open = 0
        l_atom_close = 0
        idx_string = 0
        havoc_counter = 0
        atom_index = 0
        if p_atom is not None:
            atom_index = p_atom.get_idx() + 1 + branch_index_offset

        while idx_string < len(s):
            if self._debug:
                print("\nITER={} s[idx]={}".format(idx_string, s[idx_string]))
            havoc_counter += 1
            if havoc_counter > 1000: break # just to make sure the fucking thing does not run havoc

            if s[idx_string] == "%":
                raise RingNumberNotSupported("Ring openings > 9 are currently not supported.")

            if s[idx_string] == ".":
                idx_string += 1
                if p_prev_atom is not None:
                    p_prev_atom = None
                continue

            if s[idx_string] == "(":
                i_bracket_counter += 1
                l_bracket_open = idx_string
                # now let's find the closing bracket
                for jdx_s in range(idx_string+1, len(s)):
                    if s[jdx_s] == "(":
                        i_bracket_counter += 1
                    if s[jdx_s] == ")":
                        i_bracket_counter -= 1
                    if i_bracket_counter == 0:
                        l_bracket_close = jdx_s
                        break

                if l_bracket_close == -1:
                    raise BranchBracketError
                else:
                    #print("parse branch start. indices:", l_bracket_open, l_bracket_close)
                    atom_index = self.parse_smiles(s[l_bracket_open+1:l_bracket_close], p_prev_atom, atom_index-1)
                    #print("parse branch stop. atom_index", atom_index)
                    idx_string = l_bracket_close + 1
                    continue

            if s[idx_string] == "[":
                l_atom_open = idx_string

                while s[idx_string] != "]":
                    idx_string += 1
                    if s[idx_string] == "[":
                        raise AtomBracketError
                else:
                    # found closing statement now get ready to
                    # parse the internal atom information
                    l_atom_close = idx_string
                    if self._debug:
                        print("  parse bracket atom")
                    _atoms, _bonds, new_atom_index = parse_atom(s[l_atom_open+1:l_atom_close], atom_index)

                    if p_prev_atom is not None:
                        if self._debug:
                            print("  add bond (bracket) {} <-> {}".format(p_prev_atom.get_idx(), atom_index))
                        order = 1
                        if p_bond is not None:
                            order = BOND_ORDERS[p_bond]
                            p_bond = None
                        _bond = Bond(p_prev_atom.get_idx(), atom_index, order=order)
                        self._mol.add_bond(_bond)

                    # finally add bonds internal to the parsed atom
                    self._mol.add_bonds(*_bonds)

                    p_prev_atom = _atoms[0]
                    self._mol.add_atoms(*_atoms)
                    idx_string += 1
                    atom_index = new_atom_index

            else: # atom _must_ be in the organic subset OR a specifier of sorts.
                if s[idx_string] in ["-", "=", "#", "~", ":"]:
                    # it is an explicit bond
                    p_bond = s[idx_string]
                    idx_string += 1
                    continue

                if s[idx_string] in ["B", "C", "N", "O", "P", "S", "F", "I"]: # Cl and Br handled just below
                    # check for Cl or Br
                    if self._debug:
                        print("  parse atom")
                    if s[idx_string:idx_string+2] in ["Br", "Cl"]:
                        _atoms, _bonds, _tmp = parse_atom(s[idx_string:idx_string+2], atom_index)
                        idx_string += 1
                    else:
                        _atoms, _bonds, _tmp = parse_atom(s[idx_string], atom_index)

                    # check for rings
                    try:
                        ring_index = int(s[idx_string+1]) -1 # 0 to 9
                    except ValueError: # not an integer
                        pass
                    except IndexError: # end of string reached
                        pass
                    else:
                        # here we test for ring-opening or closure
                        if self._ring_openings[ring_index] == -1:  # if closed
                            self._ring_openings[ring_index] = atom_index
                        else: # else open
                            other_atom_index = self._ring_openings[ring_index]
                            _bond = Bond(atom_index, other_atom_index)
                            self._mol.add_bond(_bond)
                            self._ring_openings[ring_index] = -1

                        idx_string += 1

                    if p_prev_atom is not None:
                        if self._debug:
                            print("  add bond {} <-> {}".format(p_prev_atom.get_idx(), atom_index))
                        order = 1
                        if p_bond is not None:
                            order = BOND_ORDERS[p_bond]
                            p_bond = None
                        _bond = Bond(p_prev_atom.get_idx(), atom_index, order=order)
                        self._mol.add_bond(_bond)

                    self._mol.add_atoms(*_atoms)
                    p_prev_atom = _atoms[0]
                    idx_string += 1
                    atom_index += 1
                    continue

        return atom_index


    def get_molecule(self):
        return self._mol


def parse_atom(s, atom_index=-1, debug=False):
    """ Parses an atom in a string s

        :param s: The string to parse
        :type s: str
        :param atom_index: the atom_index counter for continous parsing. Default is -1.
        :type atom_index: int

        :return: a list of atoms, a list of bonds and an updated atom_index for the next atom
        :rtype: list
    """
    if len(s) == 0:
        raise ValueError("parse_atom: argument 's' cannot have length 0.")

    if debug:
        print("  Smiles.parse_atom: '{}'".format(s))

    Z = 0
    if len(s) == 1:
        try:
            Z = LABEL2Z[s]
        except KeyError:
            raise IllegalAtomError("The atom '{}' is invalid.".format(s))
        else:
            # just return the atom
            return [Atom(Z, idx=atom_index)], [], atom_index +1

    idx_atom_end = -1 # atomic label from 0:idx_atom_end

    # find indices for hydrogens + counts
    n_hydrogens = 0
    idx_hydrogen = s.find("H")
    if idx_hydrogen > 0: # ignore atomic hydrogen (or proton)
        idx_atom_end = idx_hydrogen

        n_hydrogens = 1
        idx_hydrogen_count = idx_hydrogen + 1
        try:
            n_hydrogens = int(s[idx_hydrogen_count])
        except IndexError: # ran past the end of string
            pass
        except ValueError: # hit something other than a number
            pass



    idx_cat = s.find("+")
    idx_ani = s.find("-")
    idx_charge = max(idx_cat, idx_ani)
    charge = 0
    if idx_cat > 0:
        charge = 1
    elif idx_ani > 0:
        charge = -1


    if idx_charge > 0:
        if idx_hydrogen > 0:
            idx_atom_end = min(idx_charge, idx_hydrogen)
        else:
            idx_atom_end = idx_charge

        try:
            charge = int(s[idx_charge+1])
        except IndexError: # ran past the end of string
            pass
        except ValueError: # hit another + or -
            charge = charge * sum(count_items_exclusive(s, ["+", "-"]))

    if idx_atom_end == -1:
        idx_atom_end = len(s)


    if debug:
        print("  n_hydrogens :", n_hydrogens)
        print("  n_charge    :", charge)
        print("  base atom   : s[0:{}] = {}".format(idx_atom_end, s[0:idx_atom_end]))

    try:
        Z = LABEL2Z[s[0:idx_atom_end]]
    except KeyError:
        raise IllegalAtomError("The atom '{}' is invalid.".format(s[0:idx_atom_end]))

    atoms = [Atom(Z, idx=atom_index, fcharge=charge)]
    bonds = []
    for i in range(n_hydrogens):
        atoms.append(Atom(1, idx=atom_index+1+i))
        bonds.append(Bond(atom_index, atom_index+1+i))

    return atoms, bonds, atom_index+1+n_hydrogens

if __name__ == '__main__':
    SS = Smiles("C(=O)(N)c1cccnc1", debug=False)
    print(list(SS._mol.get_atoms()))
    #print(list(SS._mol.get_atoms()))
    print(list(SS._mol.get_bonds()))
