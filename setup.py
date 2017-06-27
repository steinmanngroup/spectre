#!/usr/bin/env python
#
import sys
from distutils.core import setup

# use the source code to get version information
from src.strings import version_str

__doc__="""SPECTRE: Computes spectra of molecules in solution

SPECTRE uses the polarizable embedding library in DALTON to compute
spectra for molecules in the condensed phase.
"""



def setup_spectre():
    doclines = __doc__.split("\n")

    setup(name="spectre",
        version=version_str,
        url = "https://github.com/cstein/spectre",
        author = "Casper Steinmann",
        author_email = "casper.steinmann@gmail.com",
        maintainer = "Casper Steinmann",
        maintainer_email = "casper.steinmann@gmail.com",
        license = "MIT",
        description = doclines[0],
        long_description = "\n".join(doclines[2:]),      
        #classifiers = filter(None, classifiers.split("\n")),
        platforms = ["Any."],
        package_dir={'spectre': 'src'},
        packages=['spectre', 'spectre.molecule'],
        scripts=['bin/spectre'],
        data_files=[
            ('', ['README.md','LICENSE']),
            ('share', [
                'share/dalton_loprop.bash', 'share/dalton_excited.bash'
            ])
        ]
  )

if __name__ == '__main__':
  setup_spectre()
