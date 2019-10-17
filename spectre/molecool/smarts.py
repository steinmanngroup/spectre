""" A SMARTS engine """

from molecool.smiles import parse_atom

class Smarts(object):
    """ The molecool SMARTS engine

    """

    def __init__(self, mol, s, **kwargs):
        """

        The SMARTS engine will, similarly to the SMILES engine
        generate a molecular graph.

        :param mol:
        :param s:
        :param kwargs:
        """
        self._debug = kwargs.get('debug', False)
        if self._debug:
            print("------------------------------")
            print(s)
            print("------------------------------")

        self._results = self.parse_smarts(mol, s, None)

    def parse_smarts(self, mol, s, p):
        l_idx_string = 0 # index in string
        havoc_counter = 0
        resulting_matches = []
        while l_idx_string < len(s):
            havoc_counter += 1
            if havoc_counter > 1000:
                break

            # organic subset
            if s[l_idx_string] in ["B", "C", "N", "O", "P", "S", "F", "I"]: # Cl and Br handled just below
                if s[l_idx_string:l_idx_string+1] in ["Br", "Cl"]:
                    _a, _b, _t = parse_atom(s[l_idx_string:l_idx_string+1])
                else:
                    _a, _b, _t = parse_atom(s[l_idx_string:l_idx_string+1])

                matches = match(mol, _a[0])
                if matches:
                    resulting_matches.extend(matches)

            l_idx_string += 1

        return tuple(resulting_matches[:])

def match(mol, atom):
    indices = []
    for _atom in mol.get_atoms():
        if _atom == atom:
            indices.append([_atom.get_idx()])
    return indices

