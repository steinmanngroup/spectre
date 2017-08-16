# SPECTRE
SPECTRE is a tool to compute optical properties of molecules in any homo- or heterogenous environment.

## Requirements

### Python
Preferably a version 2.7 during the initial development phase.

### Numpy
A version that ships with python. SPECTRE requires nothing fancy functionality-wise.

### OpenBabel with Python bindings
Version 2.3 or higher.
OpenBabel can be obtained from source on [github](https://github.com/openbabel/openbabel) or through your operating system package manager.

If you choose to install OpenBabel yourself remember to add your installation directory to the `PYTHONPATH` environment variable.

#### Ubuntu
On ubuntu the following should work

    sudo apt-get install openbabel python-openbabel

### FragIt
FragIt can be obtained from [github](https://github.com/FragIt/fragit-main).

Remember to add your installation directory to the `PYTHONPATH` environment variable.

### CalcIt
CalcIt can be obtained from [github](https://github.com/cstein/calcit).

Remember to add your installation directory to the `PYTHONPATH` environment variable.
Remember to export the installation directory to `CALCIT` so SPECTRE can find it.

### pepytools
peptyools can be obtained from [github](https://github.com/steinmanngroup/pepytools).

Remember to add your installation directory to the `PYTHONPATH` environment variable.

### LoProp for Dalton
LoProp for Dalton can be obtained from [github](https://github.com/vahtras/loprop).

Remember to export the installation directory to `LOPROP` so SPECTRE can find it.

### DALTON
Finally, you will need a licensed copy of the DALTON quantum chemistry program.
See [http://daltonprogram.org/](http://daltonprogram.org/) on how to obtain a license and the source code.

In order for SPECTRE to know where DALTON is you must export the `DALTON` environment variable that points to the folder where you compiled DALTON.
SPECTRE also needs to know where CALCIT is installed (slave.py script imports stuff from calcit)
For example

    export DALTON=/path/to/where/dalton/was/built

## Example
In the example folder type

    spectre -v -f pe.ini -c CCN --ex-n 4 ccn.pdb

which will compute a spectrum (excitation energies and oscillator strengths) of the chromophores `CCN` in the `ccn.pdb` file.
SPECTRE first computes the individual embedding potentials for each molecule in the `ccn.pdb` file.
Hereafter, TDDFT calculations are carried out for each chromohore embedded in the embedding potential from all other molecules.
An exciton-model is then constructed taking into account the screening effects of the environment.
Finally, a stick-spectrum is written to disk for later post-processing.
