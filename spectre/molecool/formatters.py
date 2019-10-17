from .util import Z2LABEL


class Formatter(object):
    """ Base class for formatting atoms and molecules """
    def __str__(self):
        raise NotImplemented


class AtomFormatter(Formatter):
    """ Base class for formatting atoms """
    def __init__(self, atom):
        Formatter.__init__(self)
        self._atom = atom


class MoleculeFormatter(Formatter):
    """ Base class for formatting molecules """
    def __init__(self, mol):
        Formatter.__init__(self)
        self._mol = mol


class XYZAtomFormatter(AtomFormatter):
    def __init__(self, atom):
        AtomFormatter.__init__(self, atom)

    def __str__(self):
        ATOM_LINE = "{0:<2s}{1[0]:20.9f}{1[1]:16.9f}{1[2]:16.9f}\n"
        atom_label = Z2LABEL[self._atom.get_nuclear_charge()]
        return ATOM_LINE.format(atom_label, self._atom.get_coordinate())


class XYZMoleculeFormatter(MoleculeFormatter):
    """ Formats a molecule according to the XYZ format """
    def __init__(self, mol):
        MoleculeFormatter.__init__(self, mol)


    def __str__(self):
        s_header = "{0:d}\n{1:s}\n".format(len(self._mol), self._mol.get_name())
        s_coords = ""
        for atom in self._mol.get_atoms():
            s_coords += str(XYZAtomFormatter(atom))
        return s_header + s_coords[:-1]
