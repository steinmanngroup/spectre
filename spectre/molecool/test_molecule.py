import atom
import bond
from molecule import Molecule, OBMolecule
from util import LABEL2Z

import pytest

has_openbabel = True
try:
    import openbabel
except ImportError:
    has_openbabel = False




def load_molecule_from_xyz(filename):
    mol = Molecule()
    with open(filename, 'r') as f:
        n = int(f.readline())
        title = f.readline().strip()
        mol.set_name(title)
        for i in range(n):
            tokens = f.readline().split()
            Z = LABEL2Z[tokens[0]]
            mol.add_atom(atom.Atom(Z, xyz=list(map(float, tokens[1:])), idx=i))

    return mol

def test_basic_molecule():
    """ simple tests for the Molecule class """
    mol = Molecule()
    assert mol.get_num_atoms() == 0
    assert mol.get_num_atoms() == len(mol)

    # testing basic stuff
    mol.set_name("test")
    assert mol.get_name() == "test"

    mol.set_charge(-5)
    assert mol.get_charge() == -5

    mol.set_multiplicity(1)
    assert mol.get_multiplicity() == 1

    # test that adding something that is not an atom errors out
    with pytest.raises(TypeError):
        mol.add_atom(1)



def test_molecule_add_atoms_and_bonds():
    mol = Molecule()
    _a1 = atom.Atom(1, xyz=[0.0, 0.0, 0.0], idx=0)
    _a2 = atom.Atom(1, xyz=[0.0, 0.0, 0.9], idx=1)
    mol.add_atoms(_a1, _a2)
    assert len(mol) == 2
    _b1 = bond.Bond(0, 1)
    mol.add_bonds(_b1)
    assert len(list(mol.get_bonds())) == 1

    # make sure no bond information is stored twice
    _b2 = bond.Bond(1, 0)
    mol.add_bond(_b2)
    assert len(list(mol.get_bonds())) == 1

    x,y,z = mol.get_center_of_mass()
    assert abs(z-0.45) < 1.0e-6

    with pytest.raises(IndexError):
        mol.get_atom(-1) # cannot have negative indices

    with pytest.raises(IndexError):
        mol.get_atom(100) # cannot choose atom out of range (0, 1)

    assert _a1 == mol.get_atom(0)
    assert _a2 == mol.get_atom(1)

def test_molecule_from_file():
    # test with a molecule
    mol = load_molecule_from_xyz('HOH.xyz')
    assert mol.get_num_atoms() == 3
    assert mol.get_num_atoms() == len(mol)
    assert len(list(mol.get_bonds())) == 2
    assert mol.get_name() == "HOH"


def test_molecule_copy():
    mol = load_molecule_from_xyz('CH4.xyz')
    mol2 = Molecule.from_molecule(mol)
    assert mol.get_num_atoms() == mol2.get_num_atoms()
    assert mol.get_name() == mol2.get_name()
    assert mol.get_charge() == mol2.get_charge()
    assert mol.get_multiplicity() == mol2.get_multiplicity()

    # get atoms
    assert mol.get_atom(0) == mol2.get_atom(0)
    assert len(list(mol.get_atoms())) == len(list(mol2.get_atoms()))
    for a1, a2 in zip(mol.get_atoms(), mol2.get_atoms()):
        assert a1 == a2

    # get bonds and angles
    assert len(list(mol.get_bonds())) == len(list(mol2.get_bonds()))
    assert len(list(mol.get_angles())) == len(list(mol2.get_angles()))

    # find_children
    ch1 = mol.find_children(mol.get_atom(0))
    ch2 = mol2.find_children(mol2.get_atom(0))
    assert len(ch1) == len(ch2)
    for a1, a2 in zip(ch1, ch2):
        assert a1 == a2

if has_openbabel:
    def test_obmolecule_add_atoms_and_bonds():
        mol = OBMolecule()
        _a1 = atom.Atom(1, xyz=[0.0, 0.0, 0.0], idx=0)
        _a2 = atom.Atom(1, xyz=[0.0, 0.0, 0.9], idx=1)
        mol.add_atoms(_a1, _a2)
        assert len(mol) == 2
        _b1 = bond.Bond(0, 1)
        mol.add_bonds(_b1)
        assert len(list(mol.get_bonds())) == 1

        # make sure no bond information is stored twice
        _b2 = bond.Bond(1, 0)
        mol.add_bond(_b2)
        assert len(list(mol.get_bonds())) == 1

        x,y,z = mol.get_center_of_mass()
        assert abs(z-0.45) < 1.0e-6

        with pytest.raises(IndexError):
            mol.get_atom(-1) # cannot have negative indices

        with pytest.raises(IndexError):
            mol.get_atom(100) # cannot choose atom out of range (0, 1)

        assert _a1 == mol.get_atom(0)
        assert _a2 == mol.get_atom(1)

    def test_basic_obmolecule():
        """ simple tests for the OBMolecule class """
        mol = OBMolecule()
        assert mol.get_num_atoms() == 0
        assert mol.get_num_atoms() == len(mol)

        # testing basic stuff
        mol.set_name("test")
        assert mol.get_name() == "test"

        mol.set_charge(-5)
        assert mol.get_charge() == -5

        mol.set_multiplicity(1)
        assert mol.get_multiplicity() == 1

        # test that adding something that is not an atom errors out
        with pytest.raises(TypeError):
            mol.add_atom(1)

    def test_molecule_copy_openbabel():
        """ tests if copy to openbabel molecule gives same result as regular molecule """
        mol = load_molecule_from_xyz('CH4.xyz')
        mol2 = OBMolecule.from_molecule(mol)
        assert mol.get_num_atoms() == mol2.get_num_atoms()
        assert mol.get_name() == mol2.get_name()
        assert mol.get_charge() == mol2.get_charge()
        assert mol.get_multiplicity() == mol2.get_multiplicity()

        # get atoms
        assert mol.get_atom(0) == mol2.get_atom(0)
        assert len(list(mol.get_atoms())) == len(list(mol2.get_atoms()))
        for a1, a2 in zip(mol.get_atoms(), mol2.get_atoms()):
            assert a1 == a2

        # get bonds and angles
        assert len(list(mol.get_bonds())) == len(list(mol2.get_bonds()))
        assert len(list(mol.get_angles())) == len(list(mol2.get_angles()))

        # find_children
        ch1 = mol.find_children(mol.get_atom(0))
        ch2 = mol2.find_children(mol2.get_atom(0))
        assert len(ch1) == len(ch2)
        for a1, a2 in zip(ch1, ch2):
            assert a1 == a2

    def test_bond_information_copy_correctly():
        mol1 = Molecule()
        mol1.add_atom(atom.Atom(6, xyz=[0.0, 0.0, 0.0], idx=0))
        mol1.add_atom(atom.Atom(6, xyz=[0.0, 0.0, 1.2], idx=1))
        mol1.add_bond(bond.Bond(0, 1, order=2))

        mol2 = OBMolecule.from_molecule(mol1)
        for iat, jat in zip(mol1.get_atoms(), mol2.get_atoms()):
            assert iat == jat

        for ibond, jbond in zip(mol1.get_bonds(), mol2.get_bonds()):
            assert ibond == jbond
            assert ibond.get_bond_order() == 2
            assert jbond.get_bond_order() == 2

    def test_openbabel_to_molecule_copy():
        _obmol = OBMoleculeFromFilenameAndFormat('HOH.xyz', file_format='xyz')
        mol = load_molecule_from_xyz('HOH.xyz')
        mol2 = OBMolecule(_obmol)

        assert mol.get_num_atoms() == mol2.get_num_atoms()
        assert mol.get_name() == mol2.get_name()
        assert mol.get_charge() == mol2.get_charge()
        assert mol.get_multiplicity() == mol2.get_multiplicity()

        # get atoms
        assert mol.get_atom(0) == mol2.get_atom(0)
        assert len(list(mol.get_atoms())) == len(list(mol2.get_atoms()))
        for a1, a2 in zip(mol.get_atoms(), mol2.get_atoms()):
            assert a1 == a2

        # get bonds and angles
        assert len(list(mol.get_bonds())) == len(list(mol2.get_bonds()))
        assert len(list(mol.get_angles())) == len(list(mol2.get_angles()))

        # find_children
        ch1 = mol.find_children(mol.get_atom(0))
        ch2 = mol2.find_children(mol2.get_atom(0))
        assert len(ch1) == len(ch2)
        for a1, a2 in zip(ch1, ch2):
            assert a1 == a2


    def OBMoleculeFromFilenameAndFormat(filename, file_format='pdb'):
        """ Loads a molecule into an OpenBabel molecule.

            Arguments:
            filename -- the file to load
            file_format -- file format to load

            Returns:
            OpenBabel OBMol instance of the molecule
        """
        obc = openbabel.OBConversion()
        obc.SetInFormat(file_format)
        mol = openbabel.OBMol()
        obc.ReadFile(mol, filename)
        return mol
