from __future__ import print_function

from collections import Counter
import pytest

import smiles
from smiles import Smiles

def test_smiles_parse_atoms():
    """ Tests the SMILES engine atom parser """
    SS = Smiles("") # start the engine without an atom
    assert len(SS._mol) == 0
    assert len(list(SS._mol.get_bonds())) == 0

    with pytest.raises(ValueError):
        smiles.parse_atom("")

    with pytest.raises(smiles.IllegalAtomError):
        smiles.parse_atom("Z")

    _atoms, _bonds, tmp = smiles.parse_atom("H")
    assert len(_atoms) == 1
    assert len(_bonds) == 0
    assert _atoms[0].get_nuclear_charge() == 1
    assert _atoms[0].get_formal_charge() == 0

    _atoms, _bonds, tmp = smiles.parse_atom("Cl")
    assert len(_atoms) == 1
    assert len(_bonds) == 0
    assert _atoms[0].get_nuclear_charge() == 17
    assert _atoms[0].get_formal_charge() == 0

    _atoms, _bonds, tmp = smiles.parse_atom("CH4", 0)
    assert len(_atoms) == 5
    assert len(_bonds) == 4

    _atoms, _bonds, tmp = smiles.parse_atom("H+")
    assert len(_atoms) == 1
    assert len(_bonds) == 0
    assert _atoms[0].get_nuclear_charge() == 1
    assert _atoms[0].get_formal_charge() == 1

    _atoms, _bonds, tmp = smiles.parse_atom("O-")
    assert _atoms[0].get_nuclear_charge() == 8
    assert _atoms[0].get_formal_charge() == -1

    _atoms, _bonds, tmp = smiles.parse_atom("Fe+2")
    assert _atoms[0].get_nuclear_charge() == 26
    assert _atoms[0].get_formal_charge() == 2

    _atoms, _bonds, tmp = smiles.parse_atom("Fe++")
    assert _atoms[0].get_nuclear_charge() == 26
    assert _atoms[0].get_formal_charge() == 2

    _atoms, _bonds, tmp = smiles.parse_atom("OH", 0)
    assert len(_atoms) == 2
    assert len(_bonds) == 1
    assert _atoms[0].get_nuclear_charge() == 8
    assert _atoms[0].get_formal_charge() == 0

    _atoms, _bonds, tmp = smiles.parse_atom("OH-", 0)
    assert len(_atoms) == 2
    assert len(_bonds) == 1
    assert _atoms[0].get_nuclear_charge() == 8
    assert _atoms[0].get_formal_charge() == -1

    with pytest.raises(smiles.IllegalAtomError):
        smiles.parse_atom("=N")

def test_atoms():
    SS = Smiles("F")
    _atoms = list(SS._mol.get_atoms())
    assert len(_atoms) == 1
    assert _atoms[0].get_nuclear_charge() == 9

    SS = Smiles("Cl")
    assert len(list(SS._mol.get_atoms())) == 1

    with pytest.raises(smiles.AtomBracketError):
        SS = Smiles("[C[")

    with pytest.raises(IndexError):
        SS = Smiles("[C")

    SS = Smiles("[Cl-]", debug = True)
    assert SS._mol.get_charge() == -1

    SS = Smiles("[Na+].[Cl-]", debug = True)
    assert SS._mol.get_charge() == 0

    SS = Smiles("C=[NH2+]")
    _atoms = list(SS._mol.get_atoms())
    _bonds = list(SS._mol.get_bonds())
    assert len(_atoms) == 4
    assert _bonds[0].get_bond_order() == 2
    assert SS._mol.get_charge() == 1


def test_bonds():
    SS1 = Smiles("CC")
    SS2 = Smiles("C-C")

    SS = Smiles("C=C")
    assert len(SS._mol) == 2
    bonds = list(SS._mol.get_bonds())
    assert len(bonds) == 1
    assert bonds[0].get_bond_order() == 2

    SS = Smiles("C#C")
    assert len(SS._mol) == 2
    bonds = list(SS._mol.get_bonds())
    assert len(bonds) == 1
    assert bonds[0].get_bond_order() == 3

    SS = Smiles("C=")
    assert len(SS._mol) == 1
    assert len(list(SS._mol.get_bonds())) == 0


def test_branches():
    SS = Smiles("C(C)C")
    assert len(SS._mol) == 3
    assert len(list(SS._mol.get_bonds())) == 2
    _bonds = list(SS._mol.get_bonds())
    assert _bonds[0].shares_atom(_bonds[1]) == 0 # the two bonds share the first atom

    SS = Smiles("C(=O)CC")
    assert len(SS._mol) == 4
    assert len(list(SS._mol.get_bonds())) == 3
    _bonds = list(SS._mol.get_bonds())
    assert _bonds[0].shares_atom(_bonds[1]) == 0 # the two first bonds share the first (zeroth index) atom
    assert _bonds[0].get_bond_order() == 2

    SS = Smiles("C(C(C))C")
    assert len(SS._mol) == 4
    assert len(list(SS._mol.get_bonds())) == 3
    _bonds = list(SS._mol.get_bonds())
    assert _bonds[0].shares_atom(_bonds[2]) == 0 # the first and last bonds share the first (zeroth index) atom

    with pytest.raises(smiles.BranchBracketError):
        SS = Smiles("(")

def test_rings():
    SS = Smiles("C1CC1")
    assert len(SS._mol) == 3
    assert len(list(SS._mol.get_bonds())) == 3

    SS1 = Smiles("C2CC2")
    assert len(SS._mol) == len(SS1._mol)
    assert len(list(SS._mol.get_bonds())) == len(list(SS1._mol.get_bonds()))

    with pytest.raises(smiles.RingClosureError):
        SS = Smiles("C1CC2")

    #
    # more ring openings in one string
    #
    SS = Smiles("C1C1C1C1")
    assert len(SS._mol) == 4
    assert len(list(SS._mol.get_bonds())) == 3

    # and compare with more normal example
    SS1 = Smiles("C1C1C2C2")
    assert len(SS._mol) == 4
    assert len(list(SS._mol.get_bonds())) == 3

    assert len(SS._mol) == len(SS1._mol)
    assert len(list(SS._mol.get_bonds())) == len(list(SS1._mol.get_bonds()))

    SS = Smiles("C1CC(=O)C1", debug=True)
    assert len(SS._mol) == 5
    _bonds = list(SS._mol.get_bonds())
    assert len(_bonds) == 5
    assert _bonds[2].get_bond_order() == 2

    # asserts that the bond between the two carbon
    # atoms is added only once
    SS = Smiles("C1C1")
    assert len(SS._mol) == 2
    assert len(list(SS._mol.get_bonds())) == 1
    #print(list(SS._mol.get_atoms()))
    #assert False


def test_multiple_chains():

    SS = Smiles("C(C)(C)")
    assert len(SS._mol) == 3
    _bonds = list(SS._mol.get_bonds())
    assert len(_bonds) == 2

    # explicit methane
    SS = Smiles("C([H])([H])([H])[H]")
    assert len(SS._mol) == 5
    _bonds = list(SS._mol.get_bonds())
    assert len(_bonds) == 4
    for idx, bond in enumerate(_bonds, start=1):
        assert idx == bond.get_nbr_atom_idx(0)

def test_branched_bracket_atoms():
    # the following two should provide exact same results
    SS1 = Smiles("C(O)C", debug=False)
    SS2 = Smiles("C([O])C", debug=False)

    assert len(SS1._mol) == 3
    assert len(SS1._mol) == len(SS2._mol)
    assert len(list(SS1._mol.get_bonds())) == 2
    assert len(list(SS1._mol.get_bonds())) == len(list(SS2._mol.get_bonds()))
    _atoms = list(SS2._mol.get_atoms())
    _bonds = list(SS2._mol.get_bonds())
    assert _atoms[1].get_nuclear_charge() == 8
    assert _atoms[1].get_formal_charge() == 0

    for iat, jat in zip(SS1._mol.get_atoms(), SS2._mol.get_atoms()):
        assert iat == jat

    for i, iat in enumerate(SS1._mol.get_atoms()):
        assert iat.get_idx() == i

    # test hydroxyl anion
    SS3 = Smiles("[OH1-]", debug=False)
    _atoms = list(SS3._mol.get_atoms())
    _bonds = list(SS3._mol.get_bonds())
    assert len(_atoms) == 2
    assert len(_bonds) == 1
    assert SS3._mol.get_charge() == -1



def test_nonbonded():
    SS1 = Smiles("C1.C1")
    SS2 = Smiles("CC")
    assert len(SS1._mol) == len(SS2._mol)
    _bonds1 = list(SS1._mol.get_bonds())
    _bonds2 = list(SS2._mol.get_bonds())
    assert len(_bonds1) == len(_bonds2)
    assert _bonds1[0] == _bonds2[0]

    SS1 = Smiles("[Na+].[O-]C1CCCCC1")
    SS2 = Smiles("C1CC([O-].[Na+])CCC1")
    assert len(SS1._mol) == len(SS2._mol)
    _atoms1 = list(SS1._mol.get_atoms())
    _atoms2 = list(SS2._mol.get_atoms())
    _bonds1 = list(SS1._mol.get_bonds())
    _bonds2 = list(SS2._mol.get_bonds())
    assert len(_bonds1) == len(_bonds2)

    # find the sodium atoms and make sure no bonds exist to them
    # by checking that ValueError is raised for ALL bonds in the
    # resulting molecule
    na1 = _atoms1[0]
    na2 = _atoms2[4]
    assert na1.get_nuclear_charge() == 11
    assert na2.get_nuclear_charge() == 11

    for _bond in _bonds1:
        with pytest.raises(ValueError):
            _bond.get_nbr_atom_idx(na1.get_idx())

    for _bond in _bonds2:
        with pytest.raises(ValueError):
            _bond.get_nbr_atom_idx(na2.get_idx())


def test_regressions():
    # in C([O-])C the last atom was wrong index 3 instead of 2
    # which was not found above. Looks like adding a charge to
    # oxygen makes it count indices
    SS = Smiles("C([O-])C", debug=False)
    assert len(SS._mol) == 3
    assert len(list(SS._mol.get_bonds())) == 2
    _atoms = list(SS._mol.get_atoms())
    _bonds = list(SS._mol.get_bonds())
    assert _atoms[1].get_nuclear_charge() == 8
    assert _atoms[1].get_formal_charge() == -1

    for i, iat in enumerate(SS._mol.get_atoms()):
        assert iat.get_idx() == i


