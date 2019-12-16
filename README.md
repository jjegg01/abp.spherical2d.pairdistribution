Python module abp.spherical2d.pairdistribution
==============================================

This folder contains supplementary software for the article

*J. Jeggle, J. Stenhammar and R. Wittkowski,
Pair-distribution function of active Brownian spheres in two spatial dimensions:
simulation results and analytic representation*

Contents
--------
* `demo.py`: Demo code for the module `abp.spherical2d.pairdistribution`.
See `python3 demo.py -h` for more information.
* `doc/`: HTML documentation for the module `abp.spherical2d.pairdistribution`.
* `abp/`: Python module for simplified access to the fit parameters of the
Fourier coefficients described in the article as well as routines for
reconstruction of the negative product of the pair-distribution function and the
derivative of the interaction potential. See below for installation
instructions.
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
