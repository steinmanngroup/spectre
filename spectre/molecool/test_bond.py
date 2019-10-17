import pytest

from bond import Bond


def test_creation():
    with pytest.raises(ValueError):
        Bond(0, 0)

    with pytest.raises(ValueError):
        Bond(-1, 0)

    with pytest.raises(ValueError):
        Bond(0, -1)


def test_bond():
    b = Bond(0, 1)
    assert b.get_bond_order() == 1

    b1 = Bond(1, 0)
    assert b == b1

    assert b.get_nbr_atom_idx(1) == 0


def test_atom_sharing():
    b1 = Bond(0, 1)
    b2 = Bond(1, 0)
    b3 = Bond(2, 1)
    b4 = Bond(3,2)
    assert b1.get_bond_order() == 1

    assert b1.shares_atom(b2) == -1 # they are equal

    # sharing a single atom
    assert b2.shares_atom(b3) == b3.shares_atom(b1)
    assert b3.shares_atom(b2) == b2.shares_atom(b3)

    # no sharing of atoms
    assert b4.shares_atom(b1) == -1


def test_type_checking():
    b = Bond(2,0)
    with pytest.raises(TypeError):
        value_long = float(2) # do not test long because it obsolete in python3
        b.get_nbr_atom_idx(value_long)
