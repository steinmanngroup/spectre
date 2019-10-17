from smiles import Smiles
from smarts import Smarts
import writers

def test_atomic_match():
    SS = Smiles("C", debug=False)
    mol = SS.get_molecule()
    SA = Smarts(mol, "C", debug=True)
    assert SA._results[0][0] == 0  # [0][0] -> [first match][first atom]

    SS = Smiles("CNC", debug=False)
    mol = SS.get_molecule()
    SA = Smarts(mol, "[N]", debug=True)
    assert len(SA._results) == 1  # check correct number of responses
    assert SA._results[0][0] == 1

    SA = Smarts(mol, "C", debug=True)
    assert len(SA._results) == 2  # check correct number of responses
    assert SA._results[0][0] == 0
    assert SA._results[1][0] == 2
    #print(SA._results)
    assert False
