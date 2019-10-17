molecool
========
[![Build Status](https://travis-ci.org/steinmanngroup/molecool.svg?branch=master)](https://travis-ci.org/steinmanngroup/molecool)

Molecool is a library to treated molecules treated as a
collection of atoms.
A basic use-case might be to load all atoms from an .xyz
file, manipulate the atoms and storing the molecules back
again.
The `Molecule` class can also be used as a container for
storing and manipulating atoms by other programs.

Atoms are stored in the `Atom` class which has many
properties associated with atoms.

The Molecool library does not come with file readers as
they are more likely to be out of date with what people
will use this library for, but an example of loading an
.xyz file might be:

```python
>>> with open('water.xyz', 'r') as f:
>>>     lines = f.readlines()
>>>
>>> molecule = Molecule()
>>> for atom_xyz in lines[2:]:
>>>     Z = molecule.util.LABEL2Z[atom_xyz[0]]
>>>     molecule.add_atom( Atom(Z, xyz=map(float, atom_xyz[1:])) )
>>>
>>> assert 3 == molecule.getNumAtoms()
```

The package also includes a wrapper around the OpenBabel
molecule class so it can be used as an alternative backend.

The Molecool library also comes with a [SMILES](http://www.daylight.com/dayhtml/doc/theory/theory.smiles.html)
engine which can generate molecular graphs.
It is the hope that the Molecool library will also have a [SMARTS](http://www.daylight.com/dayhtml/doc/theory/theory.smarts.html)
engine at some point.