Python module abp.spherical2d.pairdistribution
==============================================

[![DOI](https://zenodo.org/badge/228355368.svg)](https://zenodo.org/badge/latestdoi/228355368)
[![Funding provided by DFG WI 4170/3-1](https://img.shields.io/badge/DFG%20funded-WI%204170%2F3--1-blue)](https://www.dfg.de/foerderung/programme/einzelfoerderung/emmy_noether/index.html)

This folder contains supplementary software for the article

*J. Jeggle, J. Stenhammar and R. Wittkowski,
Pair-distribution function of active Brownian spheres in two spatial dimensions:
simulation results and analytic representation, J. Chem. Phys. 152, 194903 (2020),
[DOI: 10.1063/1.5140725](https://doi.org/10.1063/1.5140725)*

Contents
--------
* `abp/`: Python module for simplified access to the fit parameters of the
Fourier coefficients described in the article as well as routines for
reconstruction of the negative product of the pair-distribution function and the
derivative of the interaction potential. See below for installation
instructions.
* `demo.py`: Demo code for the module `abp.spherical2d.pairdistribution`.
See `python3 demo.py -h` for more information.
* `doc/`: HTML documentation for the module `abp.spherical2d.pairdistribution`.
* `LICENSE`: Licensing information for this software
* `README.md`: This file.

Installation
------------
Reasonably recent versions of the following software are required to make use of
the supplied code:
* Python 3
* NumPy and SciPy
* Matplotlib (needed for the demo script)

To install the Python module, copy or link the folder `abp` to a location in
your Python search path. You can find all locations in your search path by
running the following command:

```bash
python3 -c "import sys; print('\n'.join(sys.path))"
```

Please note that abp.spherical2d is a namespace package as described in PEP 420.
