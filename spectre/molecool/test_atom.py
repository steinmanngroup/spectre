import math
import pytest

import numpy

from atom import Atom


def test_atom_basic():
    a = Atom(1)
    assert a.get_nuclear_charge() == 1
    assert a.get_idx() == -1
    c = a.get_coordinate()
    assert c.dot(c) == 0.0

    with pytest.raises(TypeError):
        a.set_formal_charge(1.0)


def test_atom_coordinates():
    a = Atom(8, xyz=[0.0, 0.0, 2.0])
    assert a.get_nuclear_charge() == 8

    with pytest.raises(TypeError):
        a.set_idx(2.0)

    a.set_idx(2)
    assert a.get_idx() == 2

    r = numpy.array([0.0, 0.0, 1.0])
    with pytest.raises(TypeError):
        a.set_coordinate([])
    with pytest.raises(TypeError):
        a.set_coordinate(1.0)
    with pytest.raises(ValueError):
        a.set_coordinate(r[0:1])

    a.set_coordinate(r)
    c = a.get_coordinate()
    assert c.dot(c) == 1.0


def test_atom_coordination():
    a = Atom(6)
    assert a.get_coordination() == 4 
    with pytest.raises(TypeError):
        a.set_coordination(1.0)

    with pytest.raises(ValueError):
        a.set_coordination(5)

    a.set_coordination(3)
    assert a.get_coordination() == 3

def test_atom_hybridization():
    a = Atom(6)
    assert a.get_hybridization() == 0 # unassigned
    a.set_hybridization(3)
    assert a.get_hybridization() == 3 # assigned (here sp3 hybridized)

    with pytest.raises(TypeError):
        a.set_hybridization(1.0)

    with pytest.raises(ValueError):
        a.set_hybridization(0)

    assert a.get_hybridization() == 3 # check nothing has changed

def test_atom_copy():
    a1 = Atom(1, xyz=[0.0, 0.0, 1.0], idx=1, hybridization=3)
    a2 = Atom.from_atom(a1)

    assert a1 == a2
    assert a1.get_hybridization() == a2.get_hybridization()

    # we no longer test difference in position
    #a2 = Atom(1)
    #assert a1 != a2
