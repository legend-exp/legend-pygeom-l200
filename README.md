# legend-pygeom-l200

[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)](https://github.com/pre-commit/pre-commit)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

In this branch, we have implemented the calibration system geometry for Legend200. Calibration tubes, sources and absorbers are implemented separately from the strings, in the file called src/l200geom/calibration.py.

How to visualize it: Inside the repository, write "pip install ." and then "legend-pygeom-l200 --visualize --assemblies=calibration --fiber-modules segmented /tmp/output.gdml".
